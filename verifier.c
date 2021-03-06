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

struct verifier {
    struct puzzle_file *pf;
    struct solution_file *sf;

    int64_t throughput_cycles;
    int64_t throughput_outputs;
    int throughput_waste;

    int64_t throughput_cycles_without_poison;
    int64_t throughput_outputs_without_poison;

    uint64_t cycle_limit;

    uint64_t fails_on_wrong_output_mask;
    int wrong_output_index;
    struct board wrong_output_board;
    struct vector wrong_output_origin;
    struct vector wrong_output_basis_u;
    struct vector wrong_output_basis_v;

    int visual_loop_start_cycle;
    int visual_loop_end_cycle;

    int cycles;
    int area;

    const char *error;
};

static void *verifier_create_empty(void)
{
    struct verifier *v = calloc(sizeof(struct verifier), 1);
    v->cycle_limit = 100000;
    v->wrong_output_index = -1;
    v->cycles = -1;
    v->area = -1;
    return v;
}

void *verifier_create(const char *puzzle_filename, const char *solution_filename)
{
    struct verifier *v = verifier_create_empty();
    v->pf = parse_puzzle_file(puzzle_filename);
    if (!v->pf) {
        v->error = "invalid puzzle file";
        return v;
    }
    v->sf = parse_solution_file(solution_filename);
    if (!v->sf) {
        v->error = "invalid solution file";
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
        v->error = "invalid puzzle file";
        return v;
    }
    v->sf = parse_solution_byte_string((struct byte_string){ (unsigned char *)solution_bytes, solution_length });
    if (!v->sf) {
        v->error = "invalid solution file";
        return v;
    }
    return v;
}

const char *verifier_error(void *verifier)
{
    struct verifier *v = verifier;
    return v->error;
}

void verifier_error_clear(void *verifier)
{
    struct verifier *v = verifier;
    v->error = 0;
}

void verifier_destroy(void *verifier)
{
    struct verifier *v = verifier;
    free_puzzle_file(v->pf);
    free_solution_file(v->sf);
    verifier_wrong_output_clear(v);
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

static uint64_t min_output_count(struct solution *solution)
{
    uint64_t min = UINT64_MAX;
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        if (!(solution->inputs_and_outputs[i].type & SINGLE_OUTPUT))
            continue;
        uint64_t count = solution->inputs_and_outputs[i].number_of_outputs;
        if (count < min)
            min = count;
    }
    return min;
}

struct snapshot {
    struct mechanism *arms;
    struct board board;
    uint32_t board_in_range;
    uint64_t output_count;

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
    snapshot->output_count = min_output_count(solution);
}

static bool check_snapshot(struct solution *solution, struct board *board, struct snapshot *snapshot)
{
    if (board->collision)
        return false;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        struct mechanism a = snapshot->arms[i];
        struct mechanism b = solution->arms[i];
        if (a.position.u != b.position.u || a.position.v != b.position.v ||
         a.direction_u.u != b.direction_u.u || a.direction_u.v != b.direction_u.v)
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
        if ((a & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM)) != (b & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM)))
            return false;
        if (in_range)
            board_in_range++;
    }
    if (board_in_range != snapshot->board_in_range)
        return false;
    return true;
}

static void measure_throughput(struct verifier *v, int64_t *throughput_cycles, int64_t *throughput_outputs, int *throughput_waste, bool use_poison)
{
    struct solution solution = { 0 };
    struct board board = { 0 };
    *throughput_cycles = -1;
    *throughput_outputs = -1;
    if (!use_poison) {
        v->visual_loop_start_cycle = -1;
        v->visual_loop_end_cycle = -1;
    }
    // compute a bounding box for the solution
    int32_t max_u = INT32_MIN;
    int32_t min_u = INT32_MAX;
    int32_t max_v = INT32_MIN;
    int32_t min_v = INT32_MAX;
    for (uint32_t i = 0; i < v->sf->number_of_parts; ++i) {
        struct vector p = {
            v->sf->parts[i].position[1],
            v->sf->parts[i].position[0],
        };
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
        v->error = "no parts in solution";
        return;
    }
    // add padding of 100 units on each side
    max_u += 100;
    min_u -= 100;
    max_v += 100;
    min_v -= 100;

    if (!decode_solution(&solution, v->pf, v->sf, &v->error))
        return;
    initial_setup(&solution, &board, v->sf->area);
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
        // printf("%llu %u\n", board.cycle, throughputs_remaining);
        // print_board(&board);
        // normal throughput.
        if (!steady_state && board.cycle == snapshot.board.cycle * 2)
            take_snapshot(&solution, &board, &snapshot);
        else if (!steady_state && !(board.cycle % check_period) && snapshot.arms) {
            if (check_snapshot(&solution, &board, &snapshot)) {
                steady_state = true;
                if (!snapshot.done) {
                    snapshot.throughput_cycles = board.cycle - snapshot.board.cycle;
                    snapshot.throughput_outputs = min_output_count(&solution) - snapshot.output_count;
                    snapshot.done = true;
                    throughputs_remaining--;
                }
                if (!use_poison) {
                    v->visual_loop_start_cycle = snapshot.board.cycle;
                    v->visual_loop_end_cycle = board.cycle;
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
                v->error = "throughput measurement halted due to excessive area increase without infinite product satisfaction";
                goto error;
            }
            if (io->number_of_outputs != io->number_of_repetitions * io->outputs_per_repetition)
                continue;
            // the output is satisfied.
            if (steady_state && !--s->satisfactions_until_snapshot) {
                take_snapshot(&solution, &board, s);
                s->output_count = io->number_of_outputs;
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
                    s->throughput_outputs = io->number_of_outputs - s->output_count;
                    s->done = true;
                    throughputs_remaining--;
                }
            }
            uint32_t reps = io->number_of_repetitions + REPEATING_OUTPUT_REPETITIONS;
            if (!repeat_molecule(io, reps, &v->error))
                goto error;
        }
        cycle(&solution, &board);
    }
    if (board.collision) {
        v->error = board.collision_reason;
        goto error;
    }
    if (throughputs_remaining > 0) {
        v->error = "solution did not converge on a throughput";
        goto error;
    }
    *throughput_cycles = snapshot.throughput_cycles;
    *throughput_outputs = snapshot.throughput_outputs;
    if (throughput_waste)
        *throughput_waste = snapshot.throughput_waste;
    for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution.inputs_and_outputs[i];
        if (!(io->type & REPEATING_OUTPUT))
            continue;
        struct snapshot *s = &repeating_output_snapshots[i];
        if (s->throughput_cycles * *throughput_outputs > s->throughput_outputs * *throughput_cycles) {
            *throughput_cycles = s->throughput_cycles;
            *throughput_outputs = s->throughput_outputs;
        }
    }
error:
    for (uint32_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        free(repeating_output_snapshots[i].arms);
        free(repeating_output_snapshots[i].board.atoms_at_positions);
    }
    free(repeating_output_snapshots);
    free(shifted_atoms);
    free(snapshot.arms);
    free(snapshot.board.atoms_at_positions);
    free(shifted_board.atoms_at_positions);
    check_wrong_output_and_destroy(v, &solution, &board);
}

static void measure_dimension(struct board *board, int32_t u, int32_t v, int *dimension, int hex_width, int sortorder)
{
    int32_t max = INT32_MIN;
    int32_t min = INT32_MAX;
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID))
            continue;
        struct vector p = board->atoms_at_positions[i].position;
        int32_t value = u * p.u - v * p.v;
        if (value > max)
            max = value;
        if (value < min)
            min = value;
    }
    if (max < min)
        *dimension = 0;
    else if (sortorder * (max - min + hex_width) < sortorder * *dimension)
        *dimension = max - min + hex_width;
}

int verifier_evaluate_metric(void *verifier, const char *metric)
{
    const char *original_metric = metric;
    struct verifier *v = verifier;
    if (!v->sf)
        return -1;
    if (!strcmp(metric, "cycles") && v->cycles >= 0)
        return v->cycles;
    if (!strcmp(metric, "area (approximate)") && v->area >= 0)
        return v->area;
    long product_count = -1;
    if (!strncmp("product ", metric, strlen("product "))) {
        metric += strlen("product ");
        char *endptr = 0;
        product_count = strtol(metric, &endptr, 10);
        if (product_count < 0 || endptr == metric) {
            v->error = "invalid product count";
            return -1;
        }
        if (*endptr != ' ') {
            v->error = "product count must be followed by a metric";
            return -1;
        }
        metric = (const char *)(endptr + 1);
    }
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
                v->error = "solution refers to a reagent or product that doesn't exist in the puzzle";
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
        v->error = "invalid puzzle file";
        return -1;
    }
    if (!strcmp(metric, "throughput cycles")) {
        if (!v->throughput_cycles)
            measure_throughput(v, &v->throughput_cycles, &v->throughput_outputs, &v->throughput_waste, true);
        return v->throughput_cycles;
    } else if (!strcmp(metric, "throughput outputs")) {
        if (!v->throughput_cycles)
            measure_throughput(v, &v->throughput_cycles, &v->throughput_outputs, &v->throughput_waste, true);
        return v->throughput_outputs;
    } else if (!strcmp(metric, "throughput waste")) {
        if (!v->throughput_cycles)
            measure_throughput(v, &v->throughput_cycles, &v->throughput_outputs, &v->throughput_waste, true);
        return v->throughput_waste;
    } else if (!strcmp(metric, "throughput cycles (unrestricted)")) {
        if (!v->throughput_cycles_without_poison)
            measure_throughput(v, &v->throughput_cycles_without_poison, &v->throughput_outputs_without_poison, 0, false);
        return v->throughput_cycles_without_poison;
    } else if (!strcmp(metric, "throughput outputs (unrestricted)")) {
        if (!v->throughput_cycles_without_poison)
            measure_throughput(v, &v->throughput_cycles_without_poison, &v->throughput_outputs_without_poison, 0, false);
        return v->throughput_outputs_without_poison;
    } else if (!strcmp(metric, "visual loop start cycle")) {
        if (!v->visual_loop_start_cycle)
            measure_throughput(v, &v->throughput_cycles_without_poison, &v->throughput_outputs_without_poison, 0, false);
        return v->visual_loop_start_cycle;
    } else if (!strcmp(metric, "visual loop end cycle")) {
        if (!v->visual_loop_end_cycle)
            measure_throughput(v, &v->throughput_cycles_without_poison, &v->throughput_outputs_without_poison, 0, false);
        return v->visual_loop_end_cycle;
    }
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (!decode_solution(&solution, v->pf, v->sf, &v->error))
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
            v->error = "no hotkeys specified in 'instructions with hotkey' metric";
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
                v->error = "invalid instruction hotkey in 'instructions with hotkey' metric";
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
    if (product_count >= 0)
        solution.target_number_of_outputs = product_count;
    initial_setup(&solution, &board, v->sf->area);
    board.fails_on_wrong_output_mask = v->fails_on_wrong_output_mask;
    if (!strcmp(metric, "overlap")) {
        int overlap = INT_MAX;
        if (board.overlap < INT_MAX)
            overlap = board.overlap;
        destroy(&solution, &board);
        return overlap;
    }
    while (board.cycle < v->cycle_limit && !board.complete && !board.collision)
        cycle(&solution, &board);
    int value = -1;
    if (board.collision)
        v->error = board.collision_reason;
    else if (!board.complete)
        v->error = "solution did not complete within cycle limit";
    else {
        if (metric == original_metric) {
            v->cycles = board.cycle;
            v->area = used_area(&board);
        }
        if (!strcmp(metric, "cycles"))
            value = board.cycle;
        else if (!strcmp(metric, "area (approximate)"))
            value = used_area(&board);
        else if (!strcmp(metric, "height")) {
            value = INT_MAX;
            measure_dimension(&board, -1, 0, &value, 1, 1);
            measure_dimension(&board, 0, 1, &value, 1, 1);
            measure_dimension(&board, -1, 1, &value, 1, 1);
        } else if (!strcmp(metric, "width*2")) {
            value = INT_MAX;
            measure_dimension(&board, -1, 2, &value, 2, 1);
            measure_dimension(&board, -2, 1, &value, 2, 1);
            measure_dimension(&board, 1, 1, &value, 2, 1);
        } else if (!strcmp(metric, "omniheight")) {
            value = -INT_MAX;
            measure_dimension(&board, -1, 0, &value, 1, -1);
            measure_dimension(&board, 0, 1, &value, 1, -1);
            measure_dimension(&board, -1, 1, &value, 1, -1);
        } else if (!strcmp(metric, "omniwidth*2")) {
            value = -INT_MAX;
            measure_dimension(&board, -1, 2, &value, 2, -1);
            measure_dimension(&board, -2, 1, &value, 2, -1);
            measure_dimension(&board, 1, 1, &value, 2, -1);
        } else if (!strcmp(metric, "maximum absolute arm rotation"))
            value = solution.maximum_absolute_arm_rotation;
        else
            v->error = "unknown metric";
    }
    check_wrong_output_and_destroy(v, &solution, &board);
    return value;
}
