#include "decode.h"
#include "parse.h"
#include "sim.h"
#include "verifier.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "steady-state.h"

static void print_bonds(atom a, uint64_t bond_shift, const char *name)
{
    if (a & (1ULL << (bond_shift + 0)))
        printf(" %s[1,0]", name);
    if (a & (1ULL << (bond_shift + 1)))
        printf(" %s[0,1]", name);
    if (a & (1ULL << (bond_shift + 2)))
        printf(" %s[-1,1]", name);
    if (a & (1ULL << (bond_shift + 3)))
        printf(" %s[-1,0]", name);
    if (a & (1ULL << (bond_shift + 4)))
        printf(" %s[0,-1]", name);
    if (a & (1ULL << (bond_shift + 5)))
        printf(" %s[1,-1]", name);
}

static void print_atom(atom a)
{
    if (a & REMOVED) {
        printf(" empty");
        return;
    }
    if (a & SALT)
        printf(" salt");
    if (a & AIR)
        printf(" air");
    if (a & EARTH)
        printf(" earth");
    if (a & FIRE)
        printf(" fire");
    if (a & WATER)
        printf(" water");
    if (a & QUICKSILVER)
        printf(" quicksilver");
    if (a & GOLD)
        printf(" gold");
    if (a & SILVER)
        printf(" silver");
    if (a & COPPER)
        printf(" copper");
    if (a & IRON)
        printf(" iron");
    if (a & TIN)
        printf(" tin");
    if (a & LEAD)
        printf(" lead");
    if (a & VITAE)
        printf(" vitae");
    if (a & MORS)
        printf(" mors");
    if (a & REPEATING_OUTPUT_PLACEHOLDER)
        printf(" repetition-placeholder");
    if (a & QUINTESSENCE)
        printf(" quintessence");
    if (a & OVERLAPS_ATOMS)
        printf(" overlaps");
    if (a & BEING_PRODUCED)
        printf(" produced");
    if (a & BEING_DROPPED)
        printf(" dropped");
    if (a & VAN_BERLO_ATOM)
        printf(" van-berlo");
    if (a & IS_CHAIN_ATOM)
        printf(" chain");
    if (a & GRABBED)
        printf(" grabbed[%llu]", (a & GRABBED) / GRABBED_ONCE);
    if (a & CONDUIT_SHAPE)
        printf(" conduit");
    if (a & VISITED)
        printf(" visited");
    print_bonds(a, RECENT_BOND, "recent-bond");
    print_bonds(a, NORMAL_BOND, "bond");
    print_bonds(a, TRIPLEX_BOND_R, "triplex-r");
    print_bonds(a, TRIPLEX_BOND_Y, "triplex-y");
    print_bonds(a, TRIPLEX_BOND_K, "triplex-k");
}

static void print_board(struct board *board)
{
    int32_t maxu = -10000, minu = 10000;
    int32_t maxv = -10000, minv = 10000;
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        if (!(board->grid.atoms_at_positions[i].atom & VALID))
            continue;
        struct vector p = board->grid.atoms_at_positions[i].position;
        if (p.u < minu)
            minu = p.u;
        if (p.v < minv)
            minv = p.v;
        if (p.u > maxu)
            maxu = p.u;
        if (p.v > maxv)
            maxv = p.v;
    }
    if (maxu < minu || maxv < minv)
        return;
    int stride = (maxu - minu + 1);
    atom *points = calloc(sizeof(atom), stride * (maxv - minv + 1));
    for (uint32_t i = 0; i < BOARD_CAPACITY(board); ++i) {
        if (!(board->grid.atoms_at_positions[i].atom & VALID))
            continue;
        if (!(board->grid.atoms_at_positions[i].atom & REMOVED)) {
            printf("%d %d %"PRIx64, board->grid.atoms_at_positions[i].position.u, board->grid.atoms_at_positions[i].position.v, board->grid.atoms_at_positions[i].atom);
            print_atom(board->grid.atoms_at_positions[i].atom);
            printf("\n");
        }
        points[(board->grid.atoms_at_positions[i].position.u - minu) + stride * (board->grid.atoms_at_positions[i].position.v - minv)] = board->grid.atoms_at_positions[i].atom;
    }
    for (int v = maxv; v >= minv; --v) {
        for (int n = minv; n < v; ++n)
            printf(" ");
        for (int u = minu; u <= maxu; ++u) {
            atom a = points[stride * (v - minv) + (u - minu)];
            if (u == 0 && v == 0)
                printf(" *");
            else if (!a)
                printf("  ");
            else if (a & REMOVED)
                printf(" .");
            else if (a & IS_CHAIN_ATOM) {
                uint32_t chain = lookup_chain_atom(board, (struct vector){ u, v });
                if (chain != UINT32_MAX) {
                    if (board->chain_atoms[chain].flags & CHAIN_ATOM_IN_REPEATING_SEGMENT)
                        printf(" /");
                    else
                        printf(" %%");
                } else
                    printf(" ?");
            } else {
                for (int i = 1; i < 16; ++i) {
                    if (a & (1 << i)) {
                        printf(" %x", i & 0xf);
                        break;
                    }
                }
            }
        }
        printf("\n");
    }
    free(points);

#if 0
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->grid.atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector position = board->grid.atoms_at_positions[i].position;
        printf("%" PRId32 " %" PRId32 " %" PRIx64 "\n", position.u, position.v, a);
    }
#endif
}

int main(int argc, char *argv[])
{
    const char *puzzle_path = 0;
    const char *solution_path = 0;
    const char *requested_metrics[64];
    int number_of_requested_metrics = 0;
    bool awaiting_metric = false;
    for (int i = 1; i < argc; ++i) {
        if (awaiting_metric) {
            if (number_of_requested_metrics >= sizeof(requested_metrics) / sizeof(requested_metrics[0])) {
                fprintf(stderr, "too many metrics\n");
                return -1;
            }
            requested_metrics[number_of_requested_metrics++] = argv[i];
            awaiting_metric = false;
        } else if (strcmp(argv[i], "--metric") == 0)
            awaiting_metric = true;
        else if (!puzzle_path)
            puzzle_path = argv[i];
        else if (!solution_path)
            solution_path = argv[i];
        else {
            puzzle_path = 0;
            solution_path = 0;
            break;
        }
    }
    if (awaiting_metric || !puzzle_path || !solution_path) {
        fprintf(stderr, "usage: omsim [--metric <metric name>]* puzzle solution\n");
        return -1;
    }
    if (number_of_requested_metrics > 0) {
        void *verifier = verifier_create(puzzle_path, solution_path);
        if (verifier_error(verifier)) {
            puts(verifier_error(verifier));
            return -1;
        }
        for (int i = 0; i < number_of_requested_metrics; ++i) {
            double value = verifier_evaluate_approximate_metric(verifier, requested_metrics[i]);
            if (verifier_error(verifier)) {
                printf("%s: %s on cycle %d at %d %d\n", requested_metrics[i], verifier_error(verifier),
                    verifier_error_cycle(verifier), verifier_error_location_u(verifier), verifier_error_location_v(verifier));
                verifier_error_clear(verifier);
            } else
                printf("%s is %g\n", requested_metrics[i], value);
        }
        return 0;
    }
    // read input files.
    struct puzzle_file *pf = parse_puzzle_file(puzzle_path);
    if (!pf) {
        fprintf(stderr, "%s is not a valid puzzle file\n", puzzle_path);
        return -1;
    }
    struct solution_file *sf = parse_solution_file(solution_path);
    if (!sf) {
        fprintf(stderr, "%s is not a valid solution file\n", solution_path);
        return -1;
    }

    struct solution solution = { 0 };
    struct board board = { 0 };
    const char *error;
    if (!decode_solution(&solution, pf, sf, &error)) {
        fprintf(stderr, "solution file error: %s\n", error);
        return -1;
    }
    free_puzzle_file(pf);

    // set up the board.
    initial_setup(&solution, &board, sf->area);

    // run the solution.
    printf("-- %.*s\n", (int)sf->name.length, sf->name.bytes);
    while (board.cycle < 10000 && !board.complete) {
        printf("-- %"PRIu64" %u %u\n", board.cycle, BOARD_CAPACITY(&board), board.area);
        print_board(&board);
        cycle(&solution, &board);
        if (board.collision) {
            fprintf(stderr, "collision at %" PRId32 ", %" PRId32 ": %s\n",
             board.collision_location.u, board.collision_location.v,
             board.collision_reason);
            break;
        }
    }
    print_board(&board);
    printf("solution file says cycle count is: %" PRIu32 "\n", sf->cycles);
    printf("simulation says cycle count is: %" PRIu64 "\n", board.cycle);
    printf("solution file says area is: %" PRIu32 "\n", sf->area);
    printf("simulation says area is: %" PRIu32 "\n", used_area(&board));

#if 0
    struct steady_state steady_state = run_until_steady_state(&solution, &board, 100000);
    switch (steady_state.eventual_behavior) {
    case EVENTUALLY_ENTERS_STEADY_STATE:
        print_board(&board);
        printf("simulation entered a steady state (as of cycle %" PRIu64 "), outputting %" PRIu64 " times every %" PRIu64 " cycles\n", board.cycle, steady_state.number_of_outputs, steady_state.number_of_cycles);
        break;
    case EVENTUALLY_STOPS_RUNNING:
        printf("simulation stopped running at cycle: %" PRIu64 "\n", board.cycle);
        printf("due to collision at %" PRId32 ", %" PRId32 ": %s\n",
         board.collision_location.u, board.collision_location.v,
         board.collision_reason);
        print_board(&board);
        break;
    case EVENTUALLY_REACHES_CYCLE_LIMIT:
        printf("simulation reached cycle limit before entering a steady state\n");
        break;
    }
#endif

    destroy(&solution, &board);

    free_solution_file(sf);
    return 0;
}
