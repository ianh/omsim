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
        struct chain_atom *ca = &board->chain_atoms[i];
        // skip atoms that aren't in any list (this means they stopped being tracked as chain atoms).
        if (!ca->prev_in_list)
            continue;
        // ensure all chain atoms move in a pure translation each steady state period.
        if ((ca->flags & CHAIN_ATOM_ROTATION) != 0)
            return false;
        // avoid counting atoms in repeating segments twice.
        atom original = *lookup_atom_in_grid(&board->grid, ca->original_position);
        if ((original & VALID) && !(original & REMOVED) && !(original & IS_CHAIN_ATOM))
            ca->flags |= CHAIN_ATOM_IN_REPEATING_SEGMENT;
        else {
            ca->flags &= ~CHAIN_ATOM_IN_REPEATING_SEGMENT;
            number_of_atoms++;
        }
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

static uint64_t gcd(uint64_t a, uint64_t b)
{
    while (b != 0) {
        uint64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static struct vector chain_atom_direction(struct chain_atom *ca)
{
    int32_t u = ca->current_position.u - ca->original_position.u;
    int32_t v = ca->current_position.v - ca->original_position.v;
    int32_t divisor = gcd(labs(u), labs(v));
    return (struct vector){ u / divisor, v / divisor };
}

static int compare_chain_atoms_by_direction(const void *a, const void *b)
{
    struct vector a_dir = chain_atom_direction(*(struct chain_atom * const *)a);
    struct vector b_dir = chain_atom_direction(*(struct chain_atom * const *)b);
    if (a_dir.u - b_dir.u)
        return a_dir.u - b_dir.u;
    return a_dir.v - b_dir.v;
}

struct steady_state run_until_steady_state(struct solution *solution, struct board *board, uint64_t cycle_limit)
{
    struct snapshot snapshot = { 0 };
    uint64_t check_period = solution->tape_period;
    if (check_period == 0)
        check_period = 1;
    uint64_t next_snapshot_cycle = check_period * (1 + (board->cycle + check_period - 1) / check_period);
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        uint64_t period_aligned_start_cycle = check_period * ((solution->arm_tape_start_cycle[i] + check_period - 1) / check_period);
        if (period_aligned_start_cycle > next_snapshot_cycle)
            next_snapshot_cycle = period_aligned_start_cycle;
    }
    bool disable_check_until_next_snapshot = true;
    while (board->cycle < cycle_limit && !board->collision) {
        // printf("cycle %llu\n", board->cycle);
        if (!disable_check_until_next_snapshot && !(board->cycle % check_period) && check_snapshot(solution, board, &snapshot)) {
            // printf("check passed on cycle %llu\n", board->cycle);
            uint64_t repetition_period_length = board->cycle - snapshot.cycle;
            struct steady_state result = {
                .number_of_cycles = board->cycle - snapshot.cycle,
                .number_of_outputs = UINT64_MAX,
                .outputs_repeat_after_cycle = snapshot.cycle,
                .eventual_behavior = EVENTUALLY_ENTERS_STEADY_STATE,
            };
            for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
                struct mechanism a = snapshot.arms[i];
                struct mechanism b = solution->arms[i];
                result.pivot_parity = a.pivot_parity != b.pivot_parity;
                if (result.pivot_parity)
                    break;
            }
            for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
                if (!(solution->inputs_and_outputs[i].type & SINGLE_OUTPUT))
                    continue;
                uint64_t outputs = solution->inputs_and_outputs[i].number_of_outputs - snapshot.output_count[i];
                if (outputs < result.number_of_outputs)
                    result.number_of_outputs = outputs;
            }
            uint32_t number_of_active_chain_atoms = 0;
            uint32_t number_of_swinging_chain_atoms = 0;
            for (uint32_t i = 0; i < board->number_of_chain_atoms; ++i) {
                struct chain_atom ca = board->chain_atoms[i];
                if (!ca.prev_in_list)
                    continue;
                // chain atoms must move every period.
                if (ca.current_position.u == ca.original_position.u && ca.current_position.v == ca.original_position.v) {
                    atom *a = lookup_atom_in_grid(&board->grid, ca.current_position);
                    *a &= ~IS_CHAIN_ATOM;
                    move_chain_atom_to_list(board, i, 0);
                } else {
                    number_of_active_chain_atoms++;
                    if (ca.flags & CHAIN_ATOM_SWING_SEXTANTS)
                        number_of_swinging_chain_atoms++;
                }
            }
            if (number_of_active_chain_atoms == 0)
                board->area_growth_order = GROWTH_NONE;
            else if (number_of_swinging_chain_atoms == 0) {
                board->area_growth_order = GROWTH_LINEAR;
                for (uint32_t i = 0; i < board->number_of_area_directions; ++i)
                    free(board->area_directions[i].footprint_at_infinity.atoms_at_positions);
                free(board->area_directions);
                struct chain_atom **chain_atoms_by_direction = calloc(number_of_active_chain_atoms, sizeof(struct chain_atom *));
                number_of_active_chain_atoms = 0;
                for (uint32_t i = 0; i < board->number_of_chain_atoms; ++i) {
                    struct chain_atom *ca = &board->chain_atoms[i];
                    if (!ca->prev_in_list)
                        continue;
                    chain_atoms_by_direction[number_of_active_chain_atoms++] = ca;
                }
                qsort(chain_atoms_by_direction, number_of_active_chain_atoms, sizeof(struct chain_atom *), compare_chain_atoms_by_direction);
                board->number_of_area_directions = 0;
                struct vector previous_direction;
                for (uint32_t i = 0; i < number_of_active_chain_atoms; ++i) {
                    struct vector direction = chain_atom_direction(chain_atoms_by_direction[i]);
                    if (i == 0 || !vectors_equal(previous_direction, direction)) {
                        board->number_of_area_directions++;
                        previous_direction = direction;
                    }
                }
                board->area_directions = calloc(board->number_of_area_directions, sizeof(struct linear_area_direction));
                board->number_of_area_directions = 0;
                for (uint32_t i = 0; i < number_of_active_chain_atoms; ++i) {
                    struct vector direction = chain_atom_direction(chain_atoms_by_direction[i]);
                    if (i == 0 || !vectors_equal(previous_direction, direction)) {
                        board->area_directions[board->number_of_area_directions++].direction = direction;
                        previous_direction = direction;
                    }
                    uint32_t area_direction = board->number_of_area_directions - 1;
                    chain_atoms_by_direction[i]->area_direction = area_direction;
                    int32_t multiplier = 1;
                    if (direction.u != 0) {
                        multiplier = labs(chain_atoms_by_direction[i]->current_position.u - chain_atoms_by_direction[i]->original_position.u);
                        multiplier /= (int32_t)gcd(multiplier, labs(board->area_directions[area_direction].direction.u));
                    } else {
                        multiplier = labs(chain_atoms_by_direction[i]->current_position.v - chain_atoms_by_direction[i]->original_position.v);
                        multiplier /= (int32_t)gcd(multiplier, labs(board->area_directions[area_direction].direction.v));
                    }
                    board->area_directions[area_direction].direction.u *= multiplier;
                    board->area_directions[area_direction].direction.v *= multiplier;
                }
                free(chain_atoms_by_direction);
            } else
                board->area_growth_order = GROWTH_QUADRATIC;
            result.area_growth_order = board->area_growth_order;
            board->chain_mode = EXTEND_CHAIN;
            board->chain_will_become_visible = false;
            for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i)
                solution->inputs_and_outputs[i].maximum_feed_rate = 0;
            for (uint64_t i = 0; i < result.number_of_cycles && !board->collision; ++i) {
                // printf("last sim cycle %llu (%llu / %llu)\n", board->cycle, i, result.number_of_cycles);
                cycle(solution, board);
            }
            board->chain_mode = DISCOVER_CHAIN;
            if (board->collision)
                break;
            if (board->chain_will_become_visible) {
                // printf("re-entering box; disabling check until next snapshot\n");
                disable_check_until_next_snapshot = true;
                continue;
            }
            uint64_t repeating_outputs = 1;
            uint64_t repeating_periods = 0;
            for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
                struct input_output *io = &solution->inputs_and_outputs[i];
                if (!(io->type & REPEATING_OUTPUT))
                    continue;
                struct atom_at_position placeholder = io->original_atoms[io->number_of_original_atoms - 1];
                int32_t offset = placeholder.position.u - io->repetition_origin.u;
                if (offset <= 0 || placeholder.position.v != io->repetition_origin.v) {
                    repeating_outputs = 0;
                    repeating_periods = 1;
                    break;
                }
                if ((uint64_t)io->outputs_per_repetition * (uint64_t)io->maximum_feed_rate * repeating_periods < repeating_outputs * (uint64_t)offset) {
                    repeating_outputs = (uint64_t)io->outputs_per_repetition * (uint64_t)io->maximum_feed_rate;
                    repeating_periods = offset;
                }
            }
            if (repeating_outputs < result.number_of_outputs * repeating_periods) {
                uint64_t d = gcd(repeating_outputs, repeating_periods);
                repeating_outputs /= d;
                repeating_periods /= d;
                result.number_of_outputs = repeating_outputs;
                result.number_of_cycles *= repeating_periods;
                if (repeating_periods % 2 == 0)
                    result.pivot_parity = false;
            }
            if (board->area_growth_order == GROWTH_LINEAR) {
                uint64_t linear_area_growth = 0;
                uint64_t linear_area_growth_periods = 1;
                for (uint32_t i = 0; i < board->number_of_area_directions; ++i) {
                    for (uint32_t j = 0; j < GRID_CAPACITY(board->area_directions[i].footprint_at_infinity); ++j) {
                        struct atom_at_position ap = board->area_directions[i].footprint_at_infinity.atoms_at_positions[j];
                        if (!(ap.atom & VALID))
                            continue;
                        uint64_t divisions = ap.atom >> 1;
                        uint64_t d = gcd(divisions, linear_area_growth_periods);
                        linear_area_growth *= divisions / d;
                        linear_area_growth += linear_area_growth_periods / d;
                        linear_area_growth_periods *= divisions / d;
                    }
                }
                linear_area_growth *= result.number_of_cycles / repetition_period_length;
                uint64_t d = gcd(linear_area_growth, linear_area_growth_periods);
                linear_area_growth /= d;
                linear_area_growth_periods /= d;
                result.number_of_outputs *= linear_area_growth_periods;
                result.number_of_cycles *= linear_area_growth_periods;
                result.linear_area_growth = linear_area_growth;
            }
            destroy_snapshot(&snapshot);
            return result;
        }
        while (board->cycle > next_snapshot_cycle)
            next_snapshot_cycle *= 2;
        if (board->cycle == next_snapshot_cycle) {
            // printf("snapshot %llu\n", board->cycle);
            take_snapshot(solution, board, &snapshot);
            disable_check_until_next_snapshot = false;
            next_snapshot_cycle *= 2;
        }
        cycle(solution, board);
    }
    destroy_snapshot(&snapshot);
    return (struct steady_state){
        .eventual_behavior = board->collision ? EVENTUALLY_STOPS_RUNNING : EVENTUALLY_REACHES_CYCLE_LIMIT,
    };
}
