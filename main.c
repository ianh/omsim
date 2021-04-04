#include "parse.h"
#include "sim.h"
#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void print_board(struct board *board)
{
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector position = board->atoms_at_positions[i].position;
        printf("%" PRId32 " %" PRId32 " %" PRIx64 "\n", position.u, position.v, a);
    }
}

int main(int argc, char *argv[])
{
    lookups = 0;
    inserts = 0;
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
    if (!decode_solution(&solution, pf, sf))
        return -1;
    free_puzzle_file(pf);

    // set up the board.
    initial_setup(&solution, &board);

    // run the solution.
    printf("-- %.*s\n", (int)sf->name.length, sf->name.bytes);
    while (board.cycle < 20000 && !board.complete) {
        // printf("-- %llu %u %u\n", board.cycle, board.capacity, board.used);
        // print_board(&board);
        cycle(&solution, &board);
        if (board.collision) {
            fprintf(stderr, "collision at %" PRId32 ", %" PRId32 ": %s\n",
             board.collision_location.u, board.collision_location.v,
             board.collision_reason);
            break;
        }
    }
    printf("solution file says cycle count is: %" PRIu32 "\n", sf->cycles);
    printf("simulation says cycle count is: %" PRIu64 "\n", board.cycle);
    printf("lookups = %" PRIu64 " inserts = %" PRIu64 "\n", lookups, inserts);
    free_solution_file(sf);
    return 0;
}

