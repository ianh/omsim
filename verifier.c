#include "verifier.h"

#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

struct verifier {
    struct puzzle_file *pf;
    struct solution_file *sf;

    int throughput_cycles;
    int throughput_outputs;

    uint64_t cycle_limit;

    const char *error;
};

void *verifier_create(const char *puzzle_filename, const char *solution_filename)
{
    struct verifier *v = calloc(sizeof(struct verifier), 1);
    v->cycle_limit = 100000;
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
    free(v);
}

void verifier_set_cycle_limit(void *verifier, int cycle_limit)
{
    struct verifier *v = verifier;
    if (cycle_limit < 0)
        cycle_limit = 0;
    v->cycle_limit = cycle_limit;
}

static uint64_t min_output_count(struct solution *solution)
{
    uint64_t min = UINT64_MAX;
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        if (!(solution->inputs_and_outputs[i].type & OUTPUT))
            continue;
        uint64_t count = solution->inputs_and_outputs[i].number_of_outputs;
        if (count < min)
            min = count;
    }
    return min;
}

static void measure_throughput(struct verifier *v)
{
    struct solution solution = { 0 };
    struct board board = { 0 };
    v->throughput_cycles = -1;
    v->throughput_outputs = -1;

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
    solution.target_number_of_outputs = UINT64_MAX;
    struct mechanism *arm_snapshot = 0;
    uint64_t check_period = solution.tape_period;
    if (check_period == 0)
        check_period = 1;
    struct board board_snapshot = { .cycle = check_period };
    uint32_t board_snapshot_in_range = 0;
    uint64_t output_count_snapshot = 0;
    while (board.cycle < v->cycle_limit && !board.collision) {
        if (board.cycle == board_snapshot.cycle * 2) {
            arm_snapshot = realloc(arm_snapshot, sizeof(struct mechanism) * solution.number_of_arms);
            memcpy(arm_snapshot, solution.arms, sizeof(struct mechanism) * solution.number_of_arms);
            struct atom_at_position *a = board_snapshot.atoms_at_positions;
            board_snapshot = board;
            a = realloc(a, sizeof(struct atom_at_position) * board.capacity);
            memcpy(a, board.atoms_at_positions, sizeof(struct atom_at_position) * board.capacity);
            board_snapshot.atoms_at_positions = a;
            board_snapshot_in_range = 0;
            for (uint32_t i = 0; i < board.capacity; ++i) {
                atom a = board.atoms_at_positions[i].atom;
                if (!(a & VALID) || (a & REMOVED))
                    continue;
                struct vector p = board.atoms_at_positions[i].position;
                if (p.u < min_u || p.u > max_u || p.v < min_v || p.v > max_v)
                    continue;
                board_snapshot_in_range++;
            }
            output_count_snapshot = min_output_count(&solution);
        } else if (!(board.cycle % check_period) && arm_snapshot) {
            bool match = true;
            for (uint32_t i = 0; i < solution.number_of_arms && match; ++i) {
                struct mechanism a = arm_snapshot[i];
                struct mechanism b = solution.arms[i];
                if (a.position.u != b.position.u || a.position.v != b.position.v ||
                 a.direction_u.u != b.direction_u.u || a.direction_u.v != b.direction_u.v)
                    match = false;
            }
            uint32_t board_in_range = 0;
            for (uint32_t i = 0; i < board.capacity && match; ++i) {
                atom a = board.atoms_at_positions[i].atom;
                if (!(a & VALID) || (a & REMOVED))
                    continue;
                struct vector p = board.atoms_at_positions[i].position;
                if (p.u < min_u || p.u > max_u || p.v < min_v || p.v > max_v)
                    continue;
                atom b = *lookup_atom(&board_snapshot, p);
                if (!(b & VALID) || (b & REMOVED))
                    match = false;
                board_in_range++;
            }
            if (board_in_range != board_snapshot_in_range)
                match = false;
            if (match) {
                v->throughput_cycles = board.cycle - board_snapshot.cycle;
                v->throughput_outputs = min_output_count(&solution) - output_count_snapshot;
                break;
            }
        }
        cycle(&solution, &board);
    }
    free(board_snapshot.atoms_at_positions);
    if (board.collision)
        v->error = board.collision_reason;
    destroy(&solution, &board);
}

static void measure_dimension(struct board *board, int32_t u, int32_t v, int *dimension, int hex_width)
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
    else if (max - min + hex_width < *dimension)
        *dimension = max - min + hex_width;
}

int verifier_evaluate_metric(void *verifier, const char *metric)
{
    struct verifier *v = verifier;
    if (!v->sf) {
        v->error = "invalid solution file";
        return -1;
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
    if (!v->pf) {
        v->error = "invalid puzzle file";
        return -1;
    }
    if (!strcmp(metric, "throughput cycles")) {
        if (!v->throughput_cycles)
            measure_throughput(v);
        return v->throughput_cycles;
    } else if (!strcmp(metric, "throughput outputs")) {
        if (!v->throughput_outputs)
            measure_throughput(v);
        return v->throughput_outputs;
    }
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (!decode_solution(&solution, v->pf, v->sf, &v->error))
        return -1;
    if (!strcmp(metric, "instructions")) {
        int instructions = solution_instructions(&solution);
        destroy(&solution, &board);
        return instructions;
    } else if (!strcmp(metric, "number of arms")) {
        int arms = solution.number_of_arms;
        destroy(&solution, &board);
        return arms;
    }
    initial_setup(&solution, &board, v->sf->area);
    while (board.cycle < v->cycle_limit && !board.complete && !board.collision)
        cycle(&solution, &board);
    int value = -1;
    if (board.collision)
        v->error = board.collision_reason;
    else if (!board.complete)
        v->error = "solution did not complete within cycle limit";
    else if (!strcmp(metric, "cycles"))
        value = board.cycle;
    else if (!strcmp(metric, "area (approximate)"))
        value = used_area(&board);
    else if (!strcmp(metric, "height")) {
        value = INT_MAX;
        measure_dimension(&board, 0, -1, &value, 1);
        measure_dimension(&board, 1, 0, &value, 1);
        measure_dimension(&board, 1, -1, &value, 1);
    } else if (!strcmp(metric, "width*2")) {
        value = INT_MAX;
        measure_dimension(&board, 2, -1, &value, 2);
        measure_dimension(&board, 1, -2, &value, 2);
        measure_dimension(&board, 1, 1, &value, 2);
    } else
        v->error = "unknown metric";
    destroy(&solution, &board);
    return value;
}
