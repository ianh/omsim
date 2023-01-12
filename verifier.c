#include "verifier.h"

#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

const char *verifier_find_puzzle_name_in_solution_bytes(const char *solution_bytes,
 int solution_length, int *name_length)
{
    struct byte_string bytes = { (unsigned char *)solution_bytes, solution_length };
    struct solution_file *sf = parse_solution_byte_string(bytes);
    if (!sf)
        return 0;
    if (name_length)
        *name_length = (int)sf->puzzle.length;
    const char *name = (const char *)sf->puzzle.bytes;
    free_solution_file(sf);
    return name;
}

struct error {
    const char *description;
    int cycle;
    int location_u;
    int location_v;
};

struct per_cycle_measurements {
    int cycles;
    int area;
    int height;
    int width2;
    int omniheight;
    int omniwidth2;
    int executed_instructions;
    int maximum_absolute_arm_rotation;
    struct error error;
    bool valid;
};

struct throughput_measurements {
    int64_t throughput_cycles;
    int64_t throughput_outputs;
    int throughput_waste;
    int steady_state_start_cycle;
    int steady_state_end_cycle;
    struct per_cycle_measurements steady_state_start;
    struct per_cycle_measurements steady_state_end;
    int pivot_parity;
    struct error error;
    bool valid;
};

struct verifier {
    struct puzzle_file *pf;
    struct solution_file *sf;

    struct throughput_measurements throughput_measurements;
    struct throughput_measurements throughput_measurements_without_poison;

    uint64_t cycle_limit;

    uint64_t fails_on_wrong_output_mask;
    int wrong_output_index;
    struct board wrong_output_board;
    struct vector wrong_output_origin;
    struct vector wrong_output_basis_u;
    struct vector wrong_output_basis_v;

    struct per_cycle_measurements completion;
    struct per_cycle_measurements steady_state_start;
    struct per_cycle_measurements steady_state_end;

    int output_to_measure_intervals_for;
    int *output_intervals;
    int output_intervals_capacity;
    int number_of_output_intervals;
    int output_intervals_repeat_after;

    int throughput_margin;

    struct error error;
};

static void *verifier_create_empty(void)
{
    struct verifier *v = calloc(sizeof(struct verifier), 1);
    v->cycle_limit = 100000;
    v->wrong_output_index = -1;
    v->output_to_measure_intervals_for = -1;
    v->throughput_margin = 64;
    return v;
}

void *verifier_create(const char *puzzle_filename, const char *solution_filename)
{
    struct verifier *v = verifier_create_empty();
    v->pf = parse_puzzle_file(puzzle_filename);
    if (!v->pf) {
        v->error.description = "invalid puzzle file";
        return v;
    }
    v->sf = parse_solution_file(solution_filename);
    if (!v->sf) {
        v->error.description = "invalid solution file";
        return v;
    }
    return v;
}

void *verifier_create_from_bytes(const char *puzzle_bytes, int puzzle_length,
 const char *solution_bytes, int solution_length)
{
    char *puzzle_copy = malloc(puzzle_length);
    memcpy(puzzle_copy, puzzle_bytes, puzzle_length);
    char *solution_copy = malloc(solution_length);
    memcpy(solution_copy, solution_bytes, solution_length);
    struct verifier *v = verifier_create_from_bytes_without_copying(puzzle_copy,
     puzzle_length, solution_copy, solution_length);
    if (!v->pf || !v->sf) {
        free(puzzle_copy);
        free(solution_copy);
    } else {
        v->pf->owns_bytes = true;
        v->sf->owns_bytes = true;
    }
    return v;
}

void *verifier_create_from_bytes_without_copying(const char *puzzle_bytes, int puzzle_length,
 const char *solution_bytes, int solution_length)
{
    struct verifier *v = verifier_create_empty();
    v->pf = parse_puzzle_byte_string((struct byte_string){ (unsigned char *)puzzle_bytes, puzzle_length });
    if (!v->pf) {
        v->error.description = "invalid puzzle file";
        return v;
    }
    v->sf = parse_solution_byte_string((struct byte_string){ (unsigned char *)solution_bytes, solution_length });
    if (!v->sf) {
        v->error.description = "invalid solution file";
        return v;
    }
    return v;
}

const char *verifier_error(void *verifier)
{
    struct verifier *v = verifier;
    return v->error.description;
}

int verifier_error_cycle(void *verifier)
{
    struct verifier *v = verifier;
    return v->error.cycle;
}

int verifier_error_location_u(void *verifier)
{
    struct verifier *v = verifier;
    return v->error.location_u;
}

int verifier_error_location_v(void *verifier)
{
    struct verifier *v = verifier;
    return v->error.location_v;
}

void verifier_error_clear(void *verifier)
{
    struct verifier *v = verifier;
    v->error = (struct error){ 0 };
}

void verifier_destroy(void *verifier)
{
    struct verifier *v = verifier;
    free_puzzle_file(v->pf);
    free_solution_file(v->sf);
    verifier_wrong_output_clear(v);
    free(v->output_intervals);
    free(v);
}

void verifier_set_cycle_limit(void *verifier, int cycle_limit)
{
    struct verifier *v = verifier;
    if (cycle_limit < 0)
        cycle_limit = 0;
    v->cycle_limit = cycle_limit;
}

void verifier_set_fails_on_wrong_output(void *verifier, int output_index, int fails_on_wrong_output)
{
    struct verifier *v = verifier;
    v->fails_on_wrong_output_mask &= ~(1ULL << (uint64_t)output_index);
    v->fails_on_wrong_output_mask |= fails_on_wrong_output ? (1ULL << (uint64_t)output_index) : 0;
}
int verifier_wrong_output_index(void *verifier)
{
    struct verifier *v = verifier;
    return v->wrong_output_index;
}
int verifier_wrong_output_atom(void *verifier, int offset_u, int offset_v)
{
    struct verifier *v = verifier;
    if (v->wrong_output_index < 0)
        return -1;
    struct vector p = v->wrong_output_origin;
    p.u += v->wrong_output_basis_u.u * offset_u + v->wrong_output_basis_v.u * offset_v;
    p.v += v->wrong_output_basis_u.v * offset_u + v->wrong_output_basis_v.v * offset_v;
    atom a = *lookup_atom(&v->wrong_output_board, p);
    for (int i = 0; i <= 16; ++i) {
        if (a & (1ULL << i))
            return i;
    }
    return -1;
}
void verifier_wrong_output_clear(void *verifier)
{
    struct verifier *v = verifier;
    v->wrong_output_index = -1;
    destroy(0, &v->wrong_output_board);
}

static void check_wrong_output_and_destroy(struct verifier *v, struct solution *solution, struct board *board)
{
    if (board->wrong_output_index < solution->number_of_inputs_and_outputs && v->wrong_output_index < 0) {
        struct input_output *io = &solution->inputs_and_outputs[board->wrong_output_index];
        v->wrong_output_index = io->puzzle_index;
        v->wrong_output_board = *board;
        struct solution_part *part = &v->sf->parts[io->solution_index];
        v->wrong_output_origin = (struct vector){ part->position[0], part->position[1] };
        v->wrong_output_basis_u = u_offset_for_direction(part->rotation);
        v->wrong_output_basis_v = v_offset_for_direction(part->rotation);
        destroy(solution, 0);
    } else
        destroy(solution, board);
}

void verifier_set_throughput_margin(void *verifier, int margin)
{
    struct verifier *v = verifier;
    v->throughput_margin = margin;
    // reset throughput metrics so they're measured again if necessary.
    v->throughput_measurements = (struct throughput_measurements){ 0 };
    v->throughput_measurements_without_poison = (struct throughput_measurements){ 0 };
    v->steady_state_start = (struct per_cycle_measurements){ 0 };
    v->steady_state_end = (struct per_cycle_measurements){ 0 };
}

struct area_dimension {
    int32_t u;
    int32_t v;
    int32_t max;
    int32_t min;
};

static struct per_cycle_measurements measure_at_current_cycle(struct verifier *v, struct solution *solution, struct board *board, bool check_completion)
{
    struct per_cycle_measurements error_measurements = {
        .cycles = -1,
        .area = -1,
        .height = -1,
        .width2 = -1,
        .omniheight = -1,
        .omniwidth2 = -1,
        .executed_instructions = -1,
        .maximum_absolute_arm_rotation = -1,
        .valid = true,
    };
    if (board->collision) {
        error_measurements.error.description = board->collision_reason;
        error_measurements.error.cycle = (int)board->cycle;
        error_measurements.error.location_u = (int)board->collision_location.u;
        error_measurements.error.location_v = (int)board->collision_location.v;
        return error_measurements;
    } else if (!board->complete && check_completion) {
        error_measurements.error.description = "solution did not complete within cycle limit";
        return error_measurements;
    }
    struct area_dimension dimensions[] = {
        // height
        { -1, 0, INT32_MIN, INT32_MAX },
        { 0, 1, INT32_MIN, INT32_MAX },
        { -1, 1, INT32_MIN, INT32_MAX },

        // width
        { -1, 2, INT32_MIN, INT32_MAX },
        { -2, 1, INT32_MIN, INT32_MAX },
        { 1, 1, INT32_MIN, INT32_MAX },
    };
    bool has_atoms = false;
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID))
            continue;
        has_atoms = true;
        struct vector p = board->atoms_at_positions[i].position;
        for (int j = 0; j < sizeof(dimensions)/sizeof(dimensions[0]); ++j) {
            int32_t value = dimensions[j].u * p.u - dimensions[j].v * p.v;
            if (value < dimensions[j].min)
                dimensions[j].min = value;
            if (value > dimensions[j].max)
                dimensions[j].max = value;
        }
    }
    struct per_cycle_measurements m = {
        .cycles = (int)board->cycle,
        .area = used_area(board),
        .height = has_atoms ? INT_MAX : 0,
        .width2 = has_atoms ? INT_MAX : 0,
        .omniheight = 0,
        .omniwidth2 = 0,
        .executed_instructions = 0,
        .maximum_absolute_arm_rotation = solution ? solution->maximum_absolute_arm_rotation : -1,
        .valid = true,
    };
    if (has_atoms) {
        for (int i = 0; i < 3; ++i) {
            int32_t value = dimensions[i].max - dimensions[i].min + 1;
            if (value < m.height)
                m.height = value;
            if (value > m.omniheight)
                m.omniheight = value;
        }
        for (int i = 3; i < 6; ++i) {
            int32_t value = dimensions[i].max - dimensions[i].min + 2;
            if (value < m.width2)
                m.width2 = value;
            if (value > m.omniwidth2)
                m.omniwidth2 = value;
        }
    }
    if (solution) {
        for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
            if (board->cycle < solution->arm_tape_start_cycle[i])
                continue;
            for (size_t j = 0; j < solution->arm_tape_length[i] && j < board->cycle - solution->arm_tape_start_cycle[i]; ++j) {
                if (solution->arm_tape[i][j] != ' ' && solution->arm_tape[i][j] != '\0')
                    m.executed_instructions++;
            }
        }
    } else
        m.executed_instructions = -1;
    return m;
}

static int lookup_per_cycle_metric(struct per_cycle_measurements *measurements, const char *metric, struct error *error)
{
    if (!measurements->valid)
        return -1;
    if (!strcmp(metric, "cycles"))
        return measurements->cycles;
    else if (!strcmp(metric, "area (approximate)") || !strcmp(metric, "area"))
        return measurements->area;
    else if (!strcmp(metric, "height"))
        return measurements->height;
    else if (!strcmp(metric, "width*2"))
        return measurements->width2;
    else if (!strcmp(metric, "omniheight"))
        return measurements->omniheight;
    else if (!strcmp(metric, "omniwidth*2"))
        return measurements->omniwidth2;
    else if (!strcmp(metric, "executed instructions"))
        return measurements->executed_instructions;
    else if (!strcmp(metric, "maximum absolute arm rotation"))
        return measurements->maximum_absolute_arm_rotation;
    else {
        *error = (struct error){ .description = "unknown metric" };
        return -1;
    }
}

struct snapshot {
    struct mechanism *arms;
    struct board board;
    uint32_t board_in_range;
    uint64_t *output_count;

    int32_t max_u;
    int32_t min_u;
    int32_t max_v;
    int32_t min_v;

    bool done;

    int64_t throughput_outputs;
    int64_t throughput_cycles;
    int throughput_waste;

    // only used for repeating outputs.
    uint32_t satisfactions_until_snapshot;
    uint32_t next_satisfactions_until_snapshot;
    uint32_t number_of_repetitions;
    uint64_t last_satisfaction_cycle;
};

static void take_snapshot(struct solution *solution, struct board *board, struct snapshot *snapshot)
{
    snapshot->arms = realloc(snapshot->arms, sizeof(struct mechanism) * solution->number_of_arms);
    memcpy(snapshot->arms, solution->arms, sizeof(struct mechanism) * solution->number_of_arms);
    struct atom_at_position *a = snapshot->board.atoms_at_positions;
    snapshot->board = *board;
    a = realloc(a, sizeof(struct atom_at_position) * board->capacity);
    memcpy(a, board->atoms_at_positions, sizeof(struct atom_at_position) * board->capacity);
    snapshot->board.atoms_at_positions = a;
    snapshot->board_in_range = 0;
    snapshot->throughput_waste = 0;
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector p = board->atoms_at_positions[i].position;
        if (p.u < snapshot->min_u || p.u > snapshot->max_u || p.v < snapshot->min_v || p.v > snapshot->max_v) {
            if (board->uses_poison)
                board->atoms_at_positions[i].atom |= POISON;
            snapshot->throughput_waste = 1;
            continue;
        }
        snapshot->board_in_range++;
    }
    snapshot->output_count = realloc(snapshot->output_count, sizeof(uint64_t) * solution->number_of_inputs_and_outputs);
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i)
        snapshot->output_count[i] = solution->inputs_and_outputs[i].number_of_outputs;
}

static bool arm_directions_equivalent(struct mechanism *a, struct mechanism *b)
{
    if (a->direction_u.u == b->direction_u.u && a->direction_u.v == b->direction_u.v)
        return true;
    else if (a->type & SIX_ARM)
        return true;
    else if (a->type & TWO_ARM)
        return a->direction_u.u == -b->direction_u.u && a->direction_u.v == -b->direction_u.v;
    else if (a->type & THREE_ARM) {
        return ((-a->direction_u.u + a->direction_v.u) == b->direction_u.u && (-a->direction_u.v + a->direction_v.v) == b->direction_u.v) ||
         (-a->direction_v.u == b->direction_u.u && -a->direction_v.v == b->direction_u.v);
    } else
        return false;
}

static bool check_snapshot(struct solution *solution, struct board *board, struct snapshot *snapshot)
{
    if (board->collision)
        return false;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        struct mechanism a = snapshot->arms[i];
        struct mechanism b = solution->arms[i];
        if (a.position.u != b.position.u || a.position.v != b.position.v || !arm_directions_equivalent(&a, &b))
            return false;
    }
    uint32_t board_in_range = 0;
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector p = board->atoms_at_positions[i].position;
        bool in_range = true;
        if (p.u < snapshot->min_u || p.u > snapshot->max_u || p.v < snapshot->min_v || p.v > snapshot->max_v) {
            if (a & VISITED)
                in_range = false;
            else
                continue;
        }
        if (a & POISON)
            return false;
        atom b = *lookup_atom_without_checking_for_poison(&snapshot->board, p);
        if (!(b & VALID) || (b & REMOVED))
            return false;
        if ((a & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM | GRABBED)) != (b & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM | GRABBED)))
            return false;
        if (in_range)
            board_in_range++;
    }
    if (board_in_range != snapshot->board_in_range)
        return false;
    return true;
}

static void mark_output_interval_cycle(struct verifier *v, int cycle)
{
    if (v->number_of_output_intervals >= v->output_intervals_capacity) {
        if (v->output_intervals_capacity > 9999999) {
            v->error.description = "too many outputs while computing output intervals";
            return;
        }
        v->output_intervals_capacity = (v->output_intervals_capacity + 16) * 2;
        v->output_intervals = realloc(v->output_intervals, v->output_intervals_capacity * sizeof(v->output_intervals[0]));
    }
    v->output_intervals[v->number_of_output_intervals++] = cycle;
}

static struct throughput_measurements measure_throughput(struct verifier *v, bool use_poison)
{
    struct throughput_measurements m = {
        .valid = true,
        .throughput_cycles = -1,
        .throughput_outputs = -1,
        .steady_state_start_cycle = -1,
        .steady_state_end_cycle = -1,
    };
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (!decode_solution(&solution, v->pf, v->sf, &m.error.description))
        return m;
    initial_setup(&solution, &board, v->sf->area);
    // compute a bounding box for the solution
    int32_t max_u = INT32_MIN;
    int32_t min_u = INT32_MAX;
    int32_t max_v = INT32_MIN;
    int32_t min_v = INT32_MAX;
    for (uint32_t i = 0; i < board.capacity; ++i) {
        atom a = board.atoms_at_positions[i].atom;
        if (!(a & VALID))
            continue;
        struct vector p = board.atoms_at_positions[i].position;
        if (p.u < min_u)
            min_u = p.u;
        if (p.u > max_u)
            max_u = p.u;
        if (p.v < min_v)
            min_v = p.v;
        if (p.v > max_v)
            max_v = p.v;
    }
    if (max_u < min_u || max_v < min_v) {
        m.error.description = "no parts in solution";
        destroy(&solution, &board);
        return m;
    }
    max_u += v->throughput_margin;
    min_u -= v->throughput_margin;
    max_v += v->throughput_margin;
    min_v -= v->throughput_margin;
    board.fails_on_wrong_output_mask = v->fails_on_wrong_output_mask;
    board.ignore_swing_area = true;
    board.uses_poison = use_poison;
    board.poison_message = "solution behavior is too complex for throughput measurement";
    uint64_t check_period = solution.tape_period;
    if (check_period == 0)
        check_period = 1;
    struct snapshot snapshot = {
        .board = { .cycle = check_period },
        .max_u = max_u,
        .min_u = min_u,
        .max_v = max_v,
        .min_v = min_v,
        .throughput_cycles = 0,
        .throughput_outputs = 1,
        .throughput_waste = 0,
        .done = true,
    };
    for (uint32_t i = 0; i < solution.number_of_arms; ++i) {
        while (snapshot.board.cycle * 2 < solution.arm_tape_start_cycle[i])
            snapshot.board.cycle += check_period;
    }
    uint32_t throughputs_remaining = 0;
    struct snapshot *repeating_output_snapshots = calloc(sizeof(struct snapshot), solution.number_of_inputs_and_outputs);
    for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution.inputs_and_outputs[i];
        if (io->type & SINGLE_OUTPUT && snapshot.done) {
            throughputs_remaining++;
            snapshot.done = false;
        }
        if (!(io->type & REPEATING_OUTPUT))
            continue;
        repeating_output_snapshots[i] = (struct snapshot){
            .max_u = max_u,
            .min_u = min_u,
            .max_v = max_v,
            .min_v = min_v,
            .satisfactions_until_snapshot = 1,
            .next_satisfactions_until_snapshot = 2,
        };
        throughputs_remaining++;
    }
    solution.target_number_of_outputs = UINT64_MAX;
    struct board shifted_board = { 0 };
    struct atom_at_position *shifted_atoms = 0;
    bool steady_state = false;
    while (throughputs_remaining > 0 && board.cycle < v->cycle_limit && !board.collision) {
        // mark output intervals.
        if (!steady_state && v->output_to_measure_intervals_for >= 0) {
            for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
                if (!(solution.inputs_and_outputs[i].type & SINGLE_OUTPUT))
                    continue;
                if (solution.inputs_and_outputs[i].puzzle_index != v->output_to_measure_intervals_for)
                    continue;
                while (solution.inputs_and_outputs[i].number_of_outputs > v->number_of_output_intervals)
                    mark_output_interval_cycle(v, board.cycle);
            }
        }
        // printf("%llu %u\n", board.cycle, throughputs_remaining);
        // print_board(&board);
        // normal throughput.
        if (!steady_state && board.cycle == snapshot.board.cycle * 2) {
            v->output_intervals_repeat_after = v->number_of_output_intervals;
            take_snapshot(&solution, &board, &snapshot);
        } else if (!steady_state && !(board.cycle % check_period) && snapshot.arms) {
            if (check_snapshot(&solution, &board, &snapshot)) {
                steady_state = true;
                if (v->output_to_measure_intervals_for >= 0 && v->output_intervals_repeat_after < v->number_of_output_intervals) {
                    mark_output_interval_cycle(v, board.cycle + v->output_intervals[v->output_intervals_repeat_after] - snapshot.board.cycle);
                    v->output_intervals_repeat_after++;
                }
                if (!snapshot.done) {
                    snapshot.throughput_cycles = board.cycle - snapshot.board.cycle;
                    snapshot.throughput_outputs = 0;
                    bool outputs_set = false;
                    for (size_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
                        if (!(solution.inputs_and_outputs[i].type & SINGLE_OUTPUT))
                            continue;
                        uint64_t outputs = solution.inputs_and_outputs[i].number_of_outputs - snapshot.output_count[i];
                        if (!outputs_set || (int64_t)outputs < snapshot.throughput_outputs)
                            snapshot.throughput_outputs = outputs;
                        outputs_set = true;
                    }
                    snapshot.done = true;
                    throughputs_remaining--;
                }
                m.steady_state_start_cycle = snapshot.board.cycle;
                m.steady_state_end_cycle = board.cycle;
                m.steady_state_start = measure_at_current_cycle(v, 0, &snapshot.board, false);
                m.steady_state_end = measure_at_current_cycle(v, &solution, &board, false);
                for (uint32_t i = 0; i < solution.number_of_arms; ++i) {
                    struct mechanism a = snapshot.arms[i];
                    struct mechanism b = solution.arms[i];
                    m.pivot_parity = a.pivot_parity != b.pivot_parity;
                    if (m.pivot_parity)
                        break;
                }
            }
        }
        // repeating output throughput.
        for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
            struct input_output *io = &solution.inputs_and_outputs[i];
            if (!(io->type & REPEATING_OUTPUT))
                continue;
            struct snapshot *s = &repeating_output_snapshots[i];
            if (s->done)
                continue;
            if (board.used - s->board.used > 50000) {
                m.error.description = "throughput measurement halted due to excessive area increase without infinite product satisfaction";
                goto error;
            }
            if (io->number_of_outputs != io->number_of_repetitions * io->outputs_per_repetition) {
                // the output wasn't satisfied; periodically check whether we're spinning our wheels.
                if (steady_state && board.cycle > s->last_satisfaction_cycle + m.steady_state_end_cycle) {
                    m.error = (struct error){ .description = "throughput measurement halted due to lack of infinite product match" };
                    goto error;
                }
                continue;
            }
            // the output is satisfied.
            s->last_satisfaction_cycle = board.cycle;
            if (steady_state && !--s->satisfactions_until_snapshot) {
                take_snapshot(&solution, &board, s);
                s->number_of_repetitions = io->number_of_repetitions;
                s->satisfactions_until_snapshot = s->next_satisfactions_until_snapshot;
                s->next_satisfactions_until_snapshot *= 2;
            } else if (steady_state && s->board.cycle % check_period == board.cycle % check_period) {
                struct atom_at_position *a = shifted_board.atoms_at_positions;
                shifted_board = board;
                a = realloc(a, sizeof(struct atom_at_position) * board.capacity);
                memcpy(a, board.atoms_at_positions, sizeof(struct atom_at_position) * board.capacity);
                shifted_board.atoms_at_positions = a;

                // clear a gap corresponding to the number of repeating units added since the snapshot.
                struct atom_at_position placeholder = io->original_atoms[io->number_of_original_atoms - 1];
                struct vector offset = placeholder.position;
                offset.u -= io->repetition_origin.u;
                offset.v -= io->repetition_origin.v;
                for (uint32_t i = s->number_of_repetitions - 1; i < io->number_of_repetitions - 1; ++i) {
                    for (uint32_t j = 0; j < io->number_of_original_atoms - 1; ++j) {
                        struct vector p = io->original_atoms[j].position;
                        p.u += i * offset.u;
                        p.v += i * offset.v;
                        *lookup_atom_without_checking_for_poison(&shifted_board, p) |= REMOVED;
                    }
                }
                // shift the remaining atoms to fill the gap.
                shifted_atoms = realloc(shifted_atoms, sizeof(struct atom_at_position) * board.capacity);
                uint32_t number_of_shifted_atoms = 0;
                uint32_t shifted_atoms_cursor = 0;
                shifted_atoms[0].position = (struct vector){
                    .u = io->repetition_origin.u + (io->number_of_repetitions - 1) * offset.u,
                    .v = io->repetition_origin.v + (io->number_of_repetitions - 1) * offset.v,
                };
                atom *base = lookup_atom_without_checking_for_poison(&shifted_board, shifted_atoms[0].position);
                shifted_atoms[0].atom = *base;
                *base |= REMOVED;
                number_of_shifted_atoms++;
                while (shifted_atoms_cursor < number_of_shifted_atoms) {
                    struct atom_at_position a = shifted_atoms[shifted_atoms_cursor];
                    for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
                        if (!(a.atom & (BOND_LOW_BITS << bond_direction) & ~RECENT_BONDS))
                            continue;
                        struct vector p = a.position;
                        struct vector d = u_offset_for_direction(bond_direction);
                        p.u += d.u;
                        p.v += d.v;
                        atom *b = lookup_atom_without_checking_for_poison(&shifted_board, p);
                        if (!(*b & VALID) || (*b & REMOVED) || (*b & BEING_PRODUCED))
                            continue;
                        shifted_atoms[number_of_shifted_atoms++] = (struct atom_at_position){
                            .atom = *b,
                            .position = p,
                        };
                        *b |= REMOVED;
                    }
                    shifted_atoms_cursor++;
                }
                for (uint32_t i = 0; i < number_of_shifted_atoms; ++i) {
                    struct vector p = shifted_atoms[i].position;
                    p.u -= (io->number_of_repetitions - s->number_of_repetitions) * offset.u;
                    p.v -= (io->number_of_repetitions - s->number_of_repetitions) * offset.v;
                    // the VISITED flag tells check_snapshot not to skip over this atom, even if it's outside the bounding box.
                    *insert_atom(&shifted_board, p, "collision during shift") = (shifted_atoms[i].atom & ~POISON) | VISITED;
                }
                // compare the snapshot with the result of this gap-clearing and
                // shifting process.  if it matches, the machine is in a steady
                // state.
                if (check_snapshot(&solution, &shifted_board, s)) {
                    s->throughput_cycles = shifted_board.cycle - s->board.cycle;
                    s->throughput_outputs = io->number_of_outputs - s->output_count[i];
                    s->done = true;
                    throughputs_remaining--;
                }
            }
            uint32_t reps = io->number_of_repetitions + REPEATING_OUTPUT_REPETITIONS;
            if (!repeat_molecule(io, reps, &m.error.description))
                goto error;
        }
        cycle(&solution, &board);
    }
    if (board.collision) {
        v->output_intervals_repeat_after = -1;
        m.error.description = board.collision_reason;
        m.error.cycle = (int)board.cycle;
        m.error.location_u = (int)board.collision_location.u;
        m.error.location_v = (int)board.collision_location.v;
        goto error;
    }
    if (throughputs_remaining > 0) {
        v->output_intervals_repeat_after = -1;
        m.error.description = "solution did not converge on a throughput";
        goto error;
    }
    m.throughput_cycles = snapshot.throughput_cycles;
    m.throughput_outputs = snapshot.throughput_outputs;
    m.throughput_waste = snapshot.throughput_waste;
    for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution.inputs_and_outputs[i];
        if (!(io->type & REPEATING_OUTPUT))
            continue;
        struct snapshot *s = &repeating_output_snapshots[i];
        if (s->throughput_cycles * m.throughput_outputs > s->throughput_outputs * m.throughput_cycles) {
            m.throughput_cycles = s->throughput_cycles;
            m.throughput_outputs = s->throughput_outputs;
        }
    }
error:
    for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        free(repeating_output_snapshots[i].arms);
        free(repeating_output_snapshots[i].board.atoms_at_positions);
        free(repeating_output_snapshots[i].output_count);
    }
    free(repeating_output_snapshots);
    free(shifted_atoms);
    free(snapshot.arms);
    free(snapshot.board.atoms_at_positions);
    free(snapshot.output_count);
    free(shifted_board.atoms_at_positions);
    check_wrong_output_and_destroy(v, &solution, &board);
    return m;
}

static void ensure_output_intervals(struct verifier *v, int which_output)
{
    if (v->output_to_measure_intervals_for == which_output)
        return;

    v->output_to_measure_intervals_for = which_output;
    v->number_of_output_intervals = 0;
    v->output_intervals_repeat_after = -1;

    v->throughput_measurements = measure_throughput(v, true);

    // during measurement, the intervals are actually absolute cycles.  fix that
    // up here as a post-processing pass.
    int last = 0;
    for (int i = 0; i < v->number_of_output_intervals; ++i) {
        int delta = v->output_intervals[i] - last;
        last = v->output_intervals[i];
        v->output_intervals[i] = delta;
    }

    int cycle_start = v->output_intervals_repeat_after;
    int n = v->number_of_output_intervals - cycle_start;
    if (v->output_intervals_repeat_after > 0 && n > 0) {
        // eliminate extra repetitions within the repeating range.
        for (int i = 1; i <= n / 2; ++i) {
            if ((n % i) != 0)
                continue;
            bool repeating = true;
            for (int j = 0; j < n; ++j) {
                if (v->output_intervals[cycle_start + j] != v->output_intervals[cycle_start + (j % i)]) {
                    repeating = false;
                    break;
                }
            }
            if (repeating) {
                v->number_of_output_intervals = cycle_start + i;
                n = i;
                break;
            }
        }
        // eliminate extra repetitions before the repeating range.
        while (v->number_of_output_intervals > 0 && v->output_intervals[v->number_of_output_intervals - 1] == v->output_intervals[v->output_intervals_repeat_after - 1]) {
            v->output_intervals_repeat_after--;
            v->number_of_output_intervals--;
        }
    }
}

int verifier_number_of_output_intervals(void *verifier, int which_output)
{
    struct verifier *v = verifier;
    ensure_output_intervals(v, which_output);
    return v->number_of_output_intervals;
}

int verifier_output_interval(void *verifier, int which_output, int which_interval)
{
    struct verifier *v = verifier;
    ensure_output_intervals(v, which_output);
    if (which_interval < 0 || which_interval >= v->number_of_output_intervals)
        return -1;
    return v->output_intervals[which_interval];
}

int verifier_output_intervals_repeat_after(void *verifier, int which_output)
{
    struct verifier *v = verifier;
    ensure_output_intervals(v, which_output);
    return v->output_intervals_repeat_after;
}

int verifier_evaluate_metric(void *verifier, const char *metric)
{
    struct verifier *v = verifier;
    if (!v->sf)
        return -1;
    if (!strcmp(metric, "parsed cycles"))
        return v->sf->cycles;
    else if (!strcmp(metric, "parsed cost"))
        return v->sf->cost;
    else if (!strcmp(metric, "parsed area"))
        return v->sf->area;
    else if (!strcmp(metric, "parsed instructions"))
        return v->sf->instructions;
    else if (!strcmp(metric, "number of track segments")) {
        int value = 0;
        for (uint32_t i = 0; i < v->sf->number_of_parts; ++i)
            value += v->sf->parts[i].number_of_track_hexes;
        return value;
    } else if (!strncmp("parts of type ", metric, strlen("parts of type "))) {
        int value = 0;
        const char *part_name = metric + strlen("parts of type ");
        for (uint32_t i = 0; i < v->sf->number_of_parts; ++i) {
            if (byte_string_is(v->sf->parts[i].name, part_name))
                value++;
        }
        return value;
    } else if (!strcmp(metric, "cost"))
        return solution_file_cost(v->sf);
    else if (!strcmp(metric, "duplicate reagents") || !strcmp(metric, "duplicate products")) {
        bool input = !strcmp(metric, "duplicate reagents");
        uint32_t n = input ? v->pf->number_of_inputs : v->pf->number_of_outputs;
        int duplicates = 0;
        int *seen = calloc(n, sizeof(int));
        for (uint32_t i = 0; i < v->sf->number_of_parts; ++i) {
            if (input && !byte_string_is(v->sf->parts[i].name, "input"))
                continue;
            if (!input && !byte_string_is(v->sf->parts[i].name, "out-std") && !byte_string_is(v->sf->parts[i].name, "out-rep"))
                continue;
            uint32_t index = v->sf->parts[i].which_input_or_output;
            if (index >= n) {
                duplicates = -1;
                v->error.description = "solution refers to a reagent or product that doesn't exist in the puzzle";
                break;
            }
            if (seen[index])
                duplicates++;
            seen[index] = 1;
        }
        free(seen);
        return duplicates;
    } else if (!strcmp(metric, "maximum track gap^2")) {
        int64_t gap2 = 0;
        for (uint32_t i = 0; i < v->sf->number_of_parts; ++i) {
            for (uint32_t j = 1; j < v->sf->parts[i].number_of_track_hexes; ++j) {
                struct solution_hex_offset a = v->sf->parts[i].track_hexes[j - 1];
                struct solution_hex_offset b = v->sf->parts[i].track_hexes[j];
                int32_t du = b.offset[0] - a.offset[0];
                int32_t dv = b.offset[1] - a.offset[1];
                int64_t g = du * du + du * dv + dv * dv;
                if (g > gap2)
                    gap2 = g;
            }
        }
        return gap2 > INT_MAX ? INT_MAX : (int)gap2;
    }
    if (!v->pf) {
        v->error.description = "invalid puzzle file";
        return -1;
    }
    if (!strcmp(metric, "throughput cycles")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v, true);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_cycles;
    } else if (!strcmp(metric, "throughput outputs")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v, true);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_outputs;
    } else if (!strcmp(metric, "throughput waste")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v, true);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_waste;
    } else if (!strcmp(metric, "throughput cycles (unrestricted)")) {
        if (!v->throughput_measurements_without_poison.valid)
            v->throughput_measurements_without_poison = measure_throughput(v, false);
        v->error = v->throughput_measurements_without_poison.error;
        return v->throughput_measurements_without_poison.throughput_cycles;
    } else if (!strcmp(metric, "throughput outputs (unrestricted)")) {
        if (!v->throughput_measurements_without_poison.valid)
            v->throughput_measurements_without_poison = measure_throughput(v, false);
        v->error = v->throughput_measurements_without_poison.error;
        return v->throughput_measurements_without_poison.throughput_outputs;
    } else if (!strcmp(metric, "visual loop start cycle")) {
        if (!v->throughput_measurements_without_poison.valid)
            v->throughput_measurements_without_poison = measure_throughput(v, false);
        v->error = v->throughput_measurements_without_poison.error;
        return v->throughput_measurements_without_poison.steady_state_start_cycle;
    } else if (!strcmp(metric, "visual loop end cycle")) {
        if (!v->throughput_measurements_without_poison.valid)
            v->throughput_measurements_without_poison = measure_throughput(v, false);
        v->error = v->throughput_measurements_without_poison.error;
        int start = v->throughput_measurements_without_poison.steady_state_start_cycle;
        int end = v->throughput_measurements_without_poison.steady_state_end_cycle;
        if (v->throughput_measurements_without_poison.pivot_parity)
            return 2 * (end - start) + start;
        else
            return end;
    }
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (!decode_solution(&solution, v->pf, v->sf, &v->error.description))
        return -1;
    if (!strcmp(metric, "instructions")) {
        int instructions = solution_instructions(&solution);
        destroy(&solution, &board);
        return instructions;
    } else if (!strncmp("instructions with hotkey ", metric, strlen("instructions with hotkey "))) {
        metric += strlen("instructions with hotkey ");
        int value = 0;
        if (!*metric) {
            value = -1;
            v->error.description = "no hotkeys specified in 'instructions with hotkey' metric";
        }
        for (; *metric; ++metric) {
            switch (tolower(*metric)) {
            case 'a':
            case 'd':
            case 'e':
            case 'f':
            case 'g':
            case 'q':
            case 'r':
            case 's':
            case 't':
            case 'w':
                break;
            default:
                v->error.description = "invalid instruction hotkey in 'instructions with hotkey' metric";
                destroy(&solution, &board);
                return -1;
            }
            for (uint32_t i = 0; i < solution.number_of_arms; ++i) {
                for (size_t j = 0; j < solution.arm_tape_length[i]; ++j) {
                    if (solution.arm_tape[i][j] == tolower(*metric))
                        value++;
                }
            }
        }
        destroy(&solution, &board);
        return value;
    } else if (!strcmp(metric, "instruction tape period")) {
        int period = solution.tape_period;
        destroy(&solution, &board);
        return period;
    } else if (!strcmp(metric, "number of arms")) {
        int arms = solution.number_of_arms;
        destroy(&solution, &board);
        return arms;
    }
    initial_setup(&solution, &board, v->sf->area);
    board.fails_on_wrong_output_mask = v->fails_on_wrong_output_mask;
    if (!strcmp(metric, "overlap")) {
        int overlap = INT_MAX;
        if (board.overlap < INT_MAX)
            overlap = board.overlap;
        destroy(&solution, &board);
        return overlap;
    }
    if (!strncmp("product ", metric, strlen("product "))) {
        metric += strlen("product ");
        char *endptr = 0;
        long product_count = strtol(metric, &endptr, 10);
        if (product_count < 0 || endptr == metric) {
            v->error.description = "invalid product count";
            destroy(&solution, &board);
            return -1;
        }
        if (*endptr != ' ') {
            v->error.description = "product count must be followed by a metric";
            destroy(&solution, &board);
            return -1;
        }
        metric = (const char *)(endptr + 1);
        solution.target_number_of_outputs = product_count;
        while (board.cycle < v->cycle_limit && !board.complete && !board.collision)
            cycle(&solution, &board);
        struct per_cycle_measurements m = measure_at_current_cycle(v, &solution, &board, true);
        check_wrong_output_and_destroy(v, &solution, &board);
        v->error = m.error;
        return lookup_per_cycle_metric(&m, metric, &v->error);
    } else if (!strncmp("steady state ", metric, strlen("steady state "))) {
        metric += strlen("steady state ");
        if (!v->steady_state_start.valid || !v->steady_state_end.valid) {
            if (!v->throughput_measurements.valid)
                v->throughput_measurements = measure_throughput(v, true);
            if ((!strcmp(metric, "area (approximate)") || !strcmp(metric, "area")) && v->throughput_measurements.throughput_waste > 0) {
                v->error = (struct error){ .description = "metric doesn't reach a steady state" };
                destroy(&solution, &board);
                return -1;
            }
            if (v->throughput_measurements.throughput_cycles < 0) {
                v->error = v->throughput_measurements.error;
                destroy(&solution, &board);
                return -1;
            }
            if (lookup_per_cycle_metric(&v->throughput_measurements.steady_state_start, metric, &v->error) != lookup_per_cycle_metric(&v->throughput_measurements.steady_state_end, metric, &v->error)) {
                v->error = (struct error){ .description = "metric doesn't reach a steady state" };
                destroy(&solution, &board);
                return -1;
            }
            // measure one period later so all swings in the period are accounted for.
            int start_cycle = v->throughput_measurements.steady_state_end_cycle;
            int period = v->throughput_measurements.steady_state_end_cycle - v->throughput_measurements.steady_state_start_cycle;
            while (board.cycle < start_cycle && !board.collision)
                cycle(&solution, &board);
            v->steady_state_start = measure_at_current_cycle(v, &solution, &board, false);
            while (board.cycle < start_cycle + period && !board.collision)
                cycle(&solution, &board);
            v->steady_state_end = measure_at_current_cycle(v, &solution, &board, false);
        }
        check_wrong_output_and_destroy(v, &solution, &board);
        if (v->steady_state_start.error.description) {
            v->error = v->steady_state_start.error;
            return -1;
        }
        if (v->steady_state_end.error.description) {
            v->error = v->steady_state_end.error;
            return -1;
        }
        int start_value = lookup_per_cycle_metric(&v->steady_state_start, metric, &v->error);
        int end_value = lookup_per_cycle_metric(&v->steady_state_end, metric, &v->error);
        if (start_value < 0 || end_value < 0)
            return -1;
        if (start_value != end_value) {
            v->error = (struct error){ .description = "metric doesn't reach a steady state" };
            return -1;
        }
        return start_value;
    } else {
        if (!v->completion.valid) {
            while (board.cycle < v->cycle_limit && !board.complete && !board.collision)
                cycle(&solution, &board);
            v->completion = measure_at_current_cycle(v, &solution, &board, true);
        }
        check_wrong_output_and_destroy(v, &solution, &board);
        v->error = v->completion.error;
        return lookup_per_cycle_metric(&v->completion, metric, &v->error);
    }
}
