#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <stdio.h>
#include <stdlib.h>

struct verifier {
    struct puzzle_file *pf;
    struct solution_file *sf;

    int throughput_cycles;
    int throughput_outputs;

    const char *error;
};

void *verifier_create(const char *puzzle_filename, const char *solution_filename)
{
    struct verifier *v = calloc(sizeof(struct verifier), 1);
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

void verifier_destroy(void *verifier)
{
    struct verifier *v = verifier;
    free_puzzle_file(v->pf);
    free_solution_file(v->sf);
    free(v);
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
    if (!decode_solution(&solution, v->pf, v->sf)) {
        v->error = "unable to decode solution";
        return;
    }
    initial_setup(&solution, &board);
    solution.target_number_of_outputs = UINT64_MAX;
    struct mechanism *arm_snapshot = 0;
    struct board board_snapshot = { .cycle = solution.tape_period };
    uint32_t board_snapshot_in_range = 0;
    uint64_t output_count_snapshot = 0;
    // rough bounding box (fixme -- make this centered on the actual glyphs/arms?)
    int32_t max_u = 200;
    int32_t min_u = -200;
    int32_t max_v = 200;
    int32_t min_v = -200;
    while (board.cycle < 100000 && !board.collision) {
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
        } else if (!(board.cycle % solution.tape_period) && arm_snapshot) {
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
    if (board.collision)
        v->error = board.collision_reason;
    destroy(&solution, &board);
}

int verifier_evaluate_metric(void *verifier, const char *metric)
{
    struct verifier *v = verifier;
    if (!v->sf || !v->pf)
        return -1;
    if (!strcmp(metric, "parsed cycles"))
        return v->sf->cycles;
    else if (!strcmp(metric, "parsed cost"))
        return v->sf->cost;
    else if (!strcmp(metric, "parsed area"))
        return v->sf->area;
    else if (!strcmp(metric, "parsed instructions"))
        return v->sf->instructions;
    else if (!strcmp(metric, "cost"))
        return solution_file_cost(v->sf);
    else if (!strcmp(metric, "throughput cycles")) {
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
    if (!decode_solution(&solution, v->pf, v->sf)) {
        v->error = "unable to decode solution";
        return -1;
    }
    if (!strcmp(metric, "instructions")) {
        int instructions = solution_instructions(&solution);
        destroy(&solution, &board);
        return instructions;
    }
    initial_setup(&solution, &board);
    while (board.cycle < 100000 && !board.complete && !board.collision)
        cycle(&solution, &board);
    int value = -1;
    if (board.collision)
        v->error = board.collision_reason;
    else if (!board.complete)
        v->error = "solution did not complete within cycle limit";
    else if (!strcmp(metric, "cycles"))
        value = board.cycle;
    else
        v->error = "unknown metric";
    destroy(&solution, &board);
    return value;
}
