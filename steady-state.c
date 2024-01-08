#include "steady-state.h"

#include <stdlib.h>

// xx
#include <stdio.h>

struct snapshot {
    struct mechanism *arms;
    struct atom_grid grid;
    uint64_t number_of_atoms;
    uint64_t cycle;
    uint64_t *output_count;
};

static void take_snapshot(struct solution *solution, struct board *board, struct snapshot *snapshot)
{
    snapshot->arms = realloc(snapshot->arms, sizeof(struct mechanism) * solution->number_of_arms);
    memcpy(snapshot->arms, solution->arms, sizeof(struct mechanism) * solution->number_of_arms);
    snapshot->grid.atoms_at_positions = realloc(snapshot->grid.atoms_at_positions, sizeof(struct atom_at_position) * BOARD_CAPACITY(board));
    memcpy(snapshot->grid.atoms_at_positions, board->grid.atoms_at_positions, sizeof(struct atom_at_position) * BOARD_CAPACITY(board));
    snapshot->grid.hash_capacity = board->grid.hash_capacity;
    snapshot->cycle = board->cycle;
    snapshot->number_of_atoms = 0;
    board->number_of_chain_atoms = 0;
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        struct atom_at_position *a = &board->grid.atoms_at_positions[i];
        if (!(a->atom & VALID) || (a->atom & REMOVED))
            continue;
        a->atom &= ~IN_REPEATING_SECTION;
        if (!position_may_be_visible_to_solution(solution, a->position)) {
            a->atom |= IS_CHAIN_ATOM;
            board->number_of_chain_atoms++;
        }
        snapshot->number_of_atoms++;
    }
    snapshot->output_count = realloc(snapshot->output_count, sizeof(uint64_t) * solution->number_of_inputs_and_outputs);
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i)
        snapshot->output_count[i] = solution->inputs_and_outputs[i].number_of_outputs;

    board->chain_atoms = realloc(board->chain_atoms, sizeof(struct chain_atom) * board->number_of_chain_atoms);
    board->chain_atom_table_size = 1;
    while (board->chain_atom_table_size < board->number_of_chain_atoms)
        board->chain_atom_table_size *= 2;
    board->chain_atom_table = realloc(board->chain_atom_table, sizeof(uint32_t) * board->chain_atom_table_size);
    memset(board->chain_atom_table, 0xFF, sizeof(uint32_t) * board->chain_atom_table_size);

    uint32_t chain_atom_index = 0;
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        struct atom_at_position *a = &board->grid.atoms_at_positions[i];
        if (!(a->atom & VALID) || (a->atom & REMOVED) || !(a->atom & IS_CHAIN_ATOM))
            continue;
        board->chain_atoms[chain_atom_index] = (struct chain_atom){
            .next_in_list = UINT32_MAX,
            .current_position = a->position,
            .original_position = a->position,
        };
        add_chain_atom_to_table(board, chain_atom_index);
        chain_atom_index++;
    }
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
    uint64_t number_of_atoms = 0;
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        atom a = board->grid.atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        // chain atoms are counted separately.
        if (a & IS_CHAIN_ATOM)
            continue;
        atom b = *lookup_atom_in_grid(&snapshot->grid, board->grid.atoms_at_positions[i].position);
        if (!(b & VALID) || (b & REMOVED))
            return false;
        if ((a & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM | GRABBED)) != (b & (NORMAL_BONDS | TRIPLEX_BONDS | ANY_ATOM | GRABBED)))
            return false;
        number_of_atoms++;
    }
    for (uint32_t i = 0; i < board->number_of_chain_atoms; ++i) {
        struct chain_atom ca = board->chain_atoms[i];
        // skip atoms that aren't in any list (this means they stopped being tracked as chain atoms).
        if (!ca.prev_in_list)
            continue;
        // ensure all chain atoms move in a pure translation each steady state period.
        if (ca.rotation != 0)
            return false;
        // avoid counting atoms in repeating segments twice.
        atom original = *lookup_atom_in_grid(&board->grid, ca.original_position);
        if (!(original & VALID) || (original & REMOVED) || (original & IS_CHAIN_ATOM))
            number_of_atoms++;
    }
    if (number_of_atoms != snapshot->number_of_atoms)
        return false;
    return true;
}

static void destroy_snapshot(struct snapshot *snapshot)
{
    free(snapshot->arms);
    free(snapshot->grid.atoms_at_positions);
    free(snapshot->output_count);
}

struct steady_state run_until_steady_state(struct solution *solution, struct board *board, uint64_t cycle_limit)
{
    struct snapshot snapshot = { 0 };
    uint64_t check_period = solution->tape_period;
    if (check_period == 0)
        check_period = 1;
    uint64_t next_snapshot_cycle = check_period;
    while (board->cycle >= next_snapshot_cycle)
        next_snapshot_cycle *= 2;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        uint64_t period_aligned_start_cycle = check_period * ((solution->arm_tape_start_cycle[i] + check_period - 1) / check_period);
        if (period_aligned_start_cycle > next_snapshot_cycle)
            next_snapshot_cycle = period_aligned_start_cycle;
    }
    bool disable_check_until_next_snapshot = true;
    while (board->cycle < cycle_limit && !board->collision) {
        // printf("cycle %llu\n", board->cycle);
        while (board->cycle > next_snapshot_cycle)
            next_snapshot_cycle *= 2;
        if (board->cycle == next_snapshot_cycle) {
            take_snapshot(solution, board, &snapshot);
            disable_check_until_next_snapshot = false;
            next_snapshot_cycle *= 2;
        } else if (!disable_check_until_next_snapshot && !(board->cycle % check_period) && check_snapshot(solution, board, &snapshot)) {
            struct steady_state result = {
                .number_of_cycles = board->cycle - snapshot.cycle,
                .number_of_outputs = UINT64_MAX,
                .eventual_behavior = EVENTUALLY_ENTERS_STEADY_STATE,
            };
            for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
                if (!(solution->inputs_and_outputs[i].type & SINGLE_OUTPUT))
                    continue;
                uint64_t outputs = solution->inputs_and_outputs[i].number_of_outputs - snapshot.output_count[i];
                if (outputs < result.number_of_outputs)
                    result.number_of_outputs = outputs;
            }
            for (uint32_t i = 0; i < board->number_of_chain_atoms; ++i) {
                struct chain_atom ca = board->chain_atoms[i];
                if (!ca.prev_in_list)
                    continue;
            }
            board->chain_mode = EXTEND_CHAIN;
            board->chain_will_become_visible = false;
            for (uint64_t i = 0; i < result.number_of_cycles && !board->collision; ++i)
                cycle(solution, board);
            board->chain_mode = DISCOVER_CHAIN;
            if (board->collision)
                break;
            if (board->chain_will_become_visible) {
                disable_check_until_next_snapshot = true;
                continue;
            }
            destroy_snapshot(&snapshot);
            return result;
        }
        cycle(solution, board);
    }
    destroy_snapshot(&snapshot);
    return (struct steady_state){
        .eventual_behavior = board->collision ? EVENTUALLY_STOPS_RUNNING : EVENTUALLY_REACHES_CYCLE_LIMIT,
    };
}
