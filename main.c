#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


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
    if (a & POISON)
        printf(" poison");
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
        if (!(board->atoms_at_positions[i].atom & VALID))
            continue;
        struct vector p = board->atoms_at_positions[i].position;
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
        if (!(board->atoms_at_positions[i].atom & VALID))
            continue;
        if (!(board->atoms_at_positions[i].atom & REMOVED)) {
            printf("%d %d %llx", board->atoms_at_positions[i].position.u, board->atoms_at_positions[i].position.v, board->atoms_at_positions[i].atom);
            print_atom(board->atoms_at_positions[i].atom);
            printf("\n");
        }
        points[(board->atoms_at_positions[i].position.u - minu) + stride * (board->atoms_at_positions[i].position.v - minv)] = board->atoms_at_positions[i].atom;
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
            else {
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
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector position = board->atoms_at_positions[i].position;
        printf("%" PRId32 " %" PRId32 " %" PRIx64 "\n", position.u, position.v, a);
    }
#endif
}

int main(int argc, char *argv[])
{
    // read input files.
    if (argc < 3) {
        fprintf(stderr, "usage: omsim [puzzle] [solution]\n");
        return -1;
    }
    struct puzzle_file *pf = parse_puzzle_file(argv[1]);
    if (!pf) {
        fprintf(stderr, "%s is not a valid puzzle file\n", argv[1]);
        return -1;
    }
    struct solution_file *sf = parse_solution_file(argv[2]);
    if (!sf) {
        fprintf(stderr, "%s is not a valid solution file\n", argv[2]);
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
        printf("-- %llu %u %u\n", board.cycle, BOARD_CAPACITY(&board), board.area);
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
    destroy(&solution, &board);

    // for (int i = 0; i < solution.number_of_glyphs; ++i) {
    //     if (solution.glyphs[i].type == EQUILIBRIUM)
    //         printf("%d %d\n", solution.glyphs[i].position.u, solution.glyphs[i].position.v);
    // }
    // for (int i = 0; i < solution.number_of_inputs_and_outputs; ++i) {
    //     if (!(solution.inputs_and_outputs[i].type & INPUT))
    //         continue;
    //     printf("%d %d\n", solution.inputs_and_outputs[i].atoms[0].position.u,
    //      solution.inputs_and_outputs[i].atoms[0].position.v);
    // }

    free_solution_file(sf);
    return 0;
}

