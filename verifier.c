#include "verifier.h"

#include "decode.h"
#include "parse.h"
#include "sim.h"
#include "steady-state.h"
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

#define INSTRUCTION_D 0
#define INSTRUCTION_A 1
#define INSTRUCTION_W 2
#define INSTRUCTION_S 3
#define INSTRUCTION_R 4
#define INSTRUCTION_F 5
#define INSTRUCTION_E 6
#define INSTRUCTION_Q 7
#define INSTRUCTION_G 8
#define INSTRUCTION_T 9
#define NUMBER_OF_INSTRUCTIONS 10
static int instruction_index(char instruction)
{
    switch (instruction) {
    case 'a':
        return INSTRUCTION_A;
    case 'd':
        return INSTRUCTION_D;
    case 'e':
        return INSTRUCTION_E;
    case 'f':
        return INSTRUCTION_F;
    case 'g':
        return INSTRUCTION_G;
    case 'q':
        return INSTRUCTION_Q;
    case 'r':
        return INSTRUCTION_R;
    case 's':
        return INSTRUCTION_S;
    case 't':
        return INSTRUCTION_T;
    case 'w':
        return INSTRUCTION_W;
    default:
        return -1;
    }
}

struct per_cycle_measurements {
    int cycles;
    int area;
    int height_0;
    int height_60;
    int height_120;
    int width2_0;
    int width2_60;
    int width2_120;
    int executed_instructions;
    int instruction_executions[NUMBER_OF_INSTRUCTIONS];
    int atom_grabs[NUMBER_OF_ATOM_TYPES];
    int maximum_absolute_arm_rotation;
    struct error error;
    bool valid;
};

struct throughput_measurements {
    int64_t throughput_cycles;
    int64_t throughput_outputs;
    int64_t throughput_linear_area;
    double throughput_quadratic_area;
    enum growth_order area_growth_order;
    int throughput_waste;
    int pivot_parity;
    int steady_state_start_cycle;
    struct per_cycle_measurements steady_state;
    struct error error;
    bool valid;
};

struct verifier {
    struct puzzle_file *pf;
    struct solution_file *sf;

    struct throughput_measurements throughput_measurements;

    uint64_t cycle_limit;

    uint64_t fails_on_wrong_output_mask;
    uint64_t fails_on_wrong_output_bonds_mask;
    int wrong_output_index;
    struct board wrong_output_board;
    struct vector wrong_output_origin;
    struct vector wrong_output_basis_u;
    struct vector wrong_output_basis_v;

    struct per_cycle_measurements completion;

    int *output_intervals;
    int number_of_output_intervals;
    int output_intervals_repeat_after;

    size_t movement_limit;

    struct error error;
};

static void *verifier_create_empty(void)
{
    struct verifier *v = calloc(sizeof(struct verifier), 1);
    v->cycle_limit = 150000;
    v->movement_limit = 5000;
    v->wrong_output_index = -1;
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
void verifier_set_fails_on_wrong_output_bonds(void *verifier, int output_index, int fails_on_wrong_output_bonds)
{
    struct verifier *v = verifier;
    v->fails_on_wrong_output_bonds_mask &= ~(1ULL << (uint64_t)output_index);
    v->fails_on_wrong_output_bonds_mask |= fails_on_wrong_output_bonds ? (1ULL << (uint64_t)output_index) : 0;
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
    // this doesn't do anything since there's no throughput margin any more.
    return;
}

struct area_dimension {
    int64_t u;
    int64_t v;
    int64_t max;
    int64_t min;
};

static struct per_cycle_measurements measure_at_current_cycle(struct verifier *v, struct solution *solution, struct board *board, bool check_completion)
{
    struct per_cycle_measurements error_measurements = {
        .cycles = -1,
        .area = -1,
        .height_0 = -1,
        .height_60 = -1,
        .height_120 = -1,
        .width2_0 = -1,
        .width2_60 = -1,
        .width2_120 = -1,
        .executed_instructions = -1,
        .instruction_executions = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
        .atom_grabs = { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 },
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
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        atom a = board->grid.atoms_at_positions[i].atom;
        if (!(a & VALID))
            continue;
        has_atoms = true;
        struct vector p = board->grid.atoms_at_positions[i].position;
        for (int j = 0; j < sizeof(dimensions)/sizeof(dimensions[0]); ++j) {
            int64_t value = dimensions[j].u * p.u - dimensions[j].v * p.v;
            if (value < dimensions[j].min)
                dimensions[j].min = value;
            if (value > dimensions[j].max)
                dimensions[j].max = value;
            if (dimensions[j].max - dimensions[j].min + 2 > INT_MAX) {
                error_measurements.error.description = "solution bounding box is too large to compute width and height";
                return error_measurements;
            }
        }
    }
    struct per_cycle_measurements m = {
        .cycles = (int)board->cycle,
        .area = used_area(board),
        .executed_instructions = 0,
        .maximum_absolute_arm_rotation = solution ? solution->maximum_absolute_arm_rotation : -1,
        .valid = true,
    };
    for (int i = 0; i < NUMBER_OF_ATOM_TYPES; ++i)
        m.atom_grabs[i] = (int)board->atom_grabs[i];
    if (has_atoms) {
        m.height_0 = dimensions[0].max - dimensions[0].min + 1;
        m.height_60 = dimensions[1].max - dimensions[1].min + 1;
        m.height_120 = dimensions[2].max - dimensions[2].min + 1;
        m.width2_0 = dimensions[3].max - dimensions[3].min + 2;
        m.width2_60 = dimensions[4].max - dimensions[4].min + 2;
        m.width2_120 = dimensions[5].max - dimensions[5].min + 2;
    }
    if (solution) {
        for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
            if (board->cycle < solution->arm_tape_start_cycle[i])
                continue;
            if (solution->arm_tape_length[i] <= 0)
                continue;
            int number_of_times_through_tape = (board->cycle - solution->arm_tape_start_cycle[i]) / solution->tape_period;
            int progress_through_tape = (board->cycle - solution->arm_tape_start_cycle[i]) % solution->tape_period;
            for (size_t j = 0; j < solution->arm_tape_length[i]; ++j) {
                if (solution->arm_tape[i][j] == ' ' || solution->arm_tape[i][j] == '\0')
                    continue;
                if (j < board->cycle - solution->arm_tape_start_cycle[i])
                    m.executed_instructions++;
                m.instruction_executions[instruction_index(solution->arm_tape[i][j])] += j < progress_through_tape ? number_of_times_through_tape + 1 : number_of_times_through_tape;
            }
        }
    } else
        m.executed_instructions = -1;
    return m;
}

static const char *name_for_atom_type(int atom_type)
{
    switch (atom_type) {
    case 0:
        return "salt";
    case 1:
        return "air";
    case 2:
        return "earth";
    case 3:
        return "fire";
    case 4:
        return "water";
    case 5:
        return "quicksilver";
    case 6:
        return "gold";
    case 7:
        return "silver";
    case 8:
        return "copper";
    case 9:
        return "iron";
    case 10:
        return "tin";
    case 11:
        return "lead";
    case 12:
        return "vitae";
    case 13:
        return "mors";
    case 14:
        return 0; // ???
    case 15:
        return "quintessence";
    default:
        return 0;
    }
}

static int lookup_per_cycle_metric(struct per_cycle_measurements *measurements, const char *metric, struct error *error)
{
    if (!measurements->valid)
        return -1;
    if (!strcmp(metric, "cycles"))
        return measurements->cycles;
    else if (!strcmp(metric, "area (approximate)") || !strcmp(metric, "area"))
        return measurements->area;
    else if (!strcmp(metric, "height at 0 degrees"))
        return measurements->height_0;
    else if (!strcmp(metric, "height at 60 degrees"))
        return measurements->height_60;
    else if (!strcmp(metric, "height at 120 degrees"))
        return measurements->height_120;
    else if (!strcmp(metric, "width*2 at 0 degrees"))
        return measurements->width2_0;
    else if (!strcmp(metric, "width*2 at 60 degrees"))
        return measurements->width2_60;
    else if (!strcmp(metric, "width*2 at 120 degrees"))
        return measurements->width2_120;
    else if (!strcmp(metric, "height")) {
        int height = measurements->height_0;
        if (height < 0 || (measurements->height_60 >= 0 && measurements->height_60 < height))
            height = measurements->height_60;
        if (height < 0 || (measurements->height_120 >= 0 && measurements->height_120 < height))
            height = measurements->height_120;
        return height;
    } else if (!strcmp(metric, "width*2")) {
        int width2 = measurements->width2_0;
        if (width2 < 0 || (measurements->width2_60 >= 0 && measurements->width2_60 < width2))
            width2 = measurements->width2_60;
        if (width2 < 0 || (measurements->width2_120 >= 0 && measurements->width2_120 < width2))
            width2 = measurements->width2_120;
        return width2;
    } else if (!strcmp(metric, "executed instructions"))
        return measurements->executed_instructions;
    else if (!strcmp(metric, "maximum absolute arm rotation"))
        return measurements->maximum_absolute_arm_rotation;
    else if (!strcmp(metric, "instruction executions")) {
        int value = 0;
        for (int i = 0; i < NUMBER_OF_INSTRUCTIONS; ++i)
            value += measurements->instruction_executions[i];
        return value;
    } else if (!strncmp("instruction executions with hotkey ", metric, strlen("instruction executions with hotkey "))) {
        metric += strlen("instruction executions with hotkey ");
        if (!*metric) {
            *error = (struct error){ .description = "no hotkeys specified in 'instruction executions with hotkey' metric" };
            return -1;
        }
        int value = 0;
        for (; *metric; ++metric) {
            int idx = instruction_index(tolower(*metric));
            if (idx >= 0)
                value += measurements->instruction_executions[idx];
            else {
                *error = (struct error){ .description = "invalid instruction hotkey in 'instruction executions with hotkey' metric" };
                return -1;
            }
        }
        return value;
    } else if (!strcmp(metric, "atoms grabbed")) {
        int value = 0;
        for (int i = 0; i < NUMBER_OF_ATOM_TYPES; ++i)
            value += measurements->atom_grabs[i];
        return value;
    } else if (!strncmp("atoms grabbed of type ", metric, strlen("atoms grabbed of type "))) {
        metric += strlen("atoms grabbed of type ");
        for (int i = 0; i < NUMBER_OF_ATOM_TYPES; ++i) {
            const char *name = name_for_atom_type(i);
            if (name && !strcmp(metric, name))
                return measurements->atom_grabs[i];
        }
        *error = (struct error){ .description = "unknown atom type" };
        return -1;
    } else {
        *error = (struct error){ .description = "unknown metric" };
        return -1;
    }
}

static struct throughput_measurements measure_throughput(struct verifier *v)
{
    struct throughput_measurements m = {
        .valid = true,
        .throughput_cycles = -1,
        .throughput_outputs = -1,
        .throughput_linear_area = -1,
        .throughput_quadratic_area = -1,
        .throughput_waste = -1,
    };
    v->output_intervals_repeat_after = -1;
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (!decode_solution(&solution, v->pf, v->sf, &m.error.description))
        return m;
    initial_setup(&solution, &board, v->sf->area);
    board.movement_limit = v->movement_limit;
    board.fails_on_wrong_output_mask = v->fails_on_wrong_output_mask;
    board.fails_on_wrong_output_bonds_mask = v->fails_on_wrong_output_bonds_mask;
    struct steady_state steady_state = run_until_steady_state(&solution, &board, v->cycle_limit);

    if (board.collision) {
        m.error.description = board.collision_reason;
        m.error.cycle = (int)board.cycle;
        m.error.location_u = (int)board.collision_location.u;
        m.error.location_v = (int)board.collision_location.v;
    } else if (steady_state.eventual_behavior == EVENTUALLY_ENTERS_STEADY_STATE) {
        m.throughput_cycles = steady_state.number_of_cycles;
        m.throughput_outputs = steady_state.number_of_outputs;
        m.throughput_linear_area = steady_state.linear_area_growth;
        m.throughput_quadratic_area = steady_state.quadratic_area_growth;
        m.area_growth_order = steady_state.area_growth_order;
        m.pivot_parity = steady_state.pivot_parity;
        m.steady_state_start_cycle = steady_state.outputs_repeat_after_cycle;
        m.throughput_waste = 0;
        m.steady_state = measure_at_current_cycle(v, &solution, &board, false);
        m.steady_state.cycles = -1;
        for (uint32_t i = 0; i < board.number_of_chain_atoms; ++i) {
            struct chain_atom ca = board.chain_atoms[i];
            if (ca.prev_in_list && !vectors_equal(ca.original_position, ca.current_position)) {
                struct vector delta = {
                    ca.current_position.u - ca.original_position.u,
                    ca.current_position.v - ca.original_position.v,
                };
                bool swings = ca.flags & CHAIN_ATOM_SWING_SEXTANTS;
                if (swings || delta.u != 0)
                    m.steady_state.height_0 = -1;
                if (swings || delta.v != 0)
                    m.steady_state.height_60 = -1;
                if (swings || delta.u + delta.v != 0)
                    m.steady_state.height_120 = -1;
                if (swings || delta.u + 2 * delta.v != 0)
                    m.steady_state.width2_0 = -1;
                if (swings || 2 * delta.u + delta.v != 0)
                    m.steady_state.width2_60 = -1;
                if (swings || delta.u - delta.v != 0)
                    m.steady_state.width2_120 = -1;
                m.steady_state.area = -1;
                m.throughput_waste = 1;
            }
        }
        m.steady_state.executed_instructions = solution_instructions(&solution);
        for (uint32_t i = 0; i < solution.number_of_arms; ++i) {
            for (size_t j = 0; j < solution.arm_tape_length[i]; ++j) {
                if (solution.arm_tape[i][j] == ' ' || solution.arm_tape[i][j] == '\0')
                    continue;
                m.steady_state.instruction_executions[instruction_index(solution.arm_tape[i][j])] = -1;
            }
        }
        // xx these are unsupported right now.
        for (int i = 0; i < NUMBER_OF_ATOM_TYPES; ++i)
            m.steady_state.atom_grabs[i] = -1;
        m.steady_state.maximum_absolute_arm_rotation = -1;
    } else
        m.error.description = "solution did not converge on a throughput";

    // output intervals are currently unsupported for puzzles with polymers.
    bool output_intervals_supported = true;
    for (size_t i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
        if (solution.inputs_and_outputs[i].type & REPEATING_OUTPUT) {
            output_intervals_supported = false;
            break;
        }
    }
    if (output_intervals_supported) {
        v->output_intervals = calloc(board.number_of_output_cycles, sizeof(uint64_t));
        for (uint64_t i = 0; i < board.number_of_output_cycles; ++i)
            v->output_intervals[i] = board.output_cycles[i];
        v->number_of_output_intervals = board.number_of_output_cycles;

        // if the solution enters a steady state, outputs repeat during the steady state period.
        if (steady_state.eventual_behavior == EVENTUALLY_ENTERS_STEADY_STATE) {
            for (size_t i = 0; i < v->number_of_output_intervals; ++i) {
                if (v->output_intervals[i] > steady_state.outputs_repeat_after_cycle) {
                    v->output_intervals_repeat_after = i;
                    break;
                }
            }
        }

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
    check_wrong_output_and_destroy(v, &solution, &board);
    return m;
}

static void ensure_output_intervals(struct verifier *v)
{
    if (v->throughput_measurements.valid)
        return;
    v->throughput_measurements = measure_throughput(v);
}

int verifier_number_of_output_intervals(void *verifier)
{
    struct verifier *v = verifier;
    if (!v->sf)
        return 0;
    ensure_output_intervals(v);
    return v->number_of_output_intervals;
}

int verifier_output_interval(void *verifier, int which_interval)
{
    struct verifier *v = verifier;
    if (!v->sf)
        return -1;
    ensure_output_intervals(v);
    if (which_interval < 0 || which_interval >= v->number_of_output_intervals)
        return -1;
    return v->output_intervals[which_interval];
}

int verifier_output_intervals_repeat_after(void *verifier)
{
    struct verifier *v = verifier;
    if (!v->sf)
        return -1;
    ensure_output_intervals(v);
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
    if (!strcmp(metric, "reaches steady state")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        if (v->throughput_measurements.throughput_cycles < 0)
            return 0;
        else
            return 1;
    } else if (!strcmp(metric, "per repetition cycles") || !strcmp(metric, "throughput cycles")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_cycles;
    } else if (!strcmp(metric, "per repetition outputs") || !strcmp(metric, "throughput outputs")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_outputs;
    } else if (!strcmp(metric, "per repetition area")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        if (v->throughput_measurements.area_growth_order == GROWTH_QUADRATIC) {
            v->error = (struct error){ .description = "metric doesn't have a consistent per-repetition value" };
            return -1;
        }
        return v->throughput_measurements.throughput_linear_area;
    } else if (!strcmp(metric, "throughput waste")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_waste;
    } else if (!strcmp(metric, "visual loop start cycle")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.steady_state_start_cycle;
    } else if (!strcmp(metric, "visual loop end cycle")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        int start = v->throughput_measurements.steady_state_start_cycle;
        int period = v->throughput_measurements.throughput_cycles;
        if (v->throughput_measurements.pivot_parity)
            return start + 2 * period;
        else
            return start + period;
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
    board.movement_limit = v->movement_limit;
    board.fails_on_wrong_output_mask = v->fails_on_wrong_output_mask;
    board.fails_on_wrong_output_bonds_mask = v->fails_on_wrong_output_bonds_mask;
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
        destroy(&solution, &board);
        metric += strlen("steady state ");
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        if (v->throughput_measurements.throughput_cycles < 0) {
            v->error = v->throughput_measurements.error;
            return -1;
        }
        v->error = v->throughput_measurements.steady_state.error;
        int value = lookup_per_cycle_metric(&v->throughput_measurements.steady_state, metric, &v->error);
        if (value < 0 && !v->error.description)
            v->error = (struct error){ .description = "metric doesn't reach a steady state" };
        return value;
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

double verifier_evaluate_approximate_metric(void *verifier, const char *metric)
{
    struct verifier *v = verifier;
    if (!v->sf)
        return -1;
    if (!strcmp(metric, "per repetition^2 area")) {
        if (!v->throughput_measurements.valid)
            v->throughput_measurements = measure_throughput(v);
        v->error = v->throughput_measurements.error;
        return v->throughput_measurements.throughput_quadratic_area;
    } else
        return verifier_evaluate_metric(verifier, metric);
}
