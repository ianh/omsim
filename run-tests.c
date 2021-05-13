#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

struct puzzle {
    struct puzzle *next;
    char *filename;
    struct puzzle_file *pf;
};

int main()
{
    size_t n = 64;
    char *buf = malloc(n);

    struct puzzle *puzzles = 0;
    FILE *puzzle_list = popen("find test/puzzle -type f", "r");
    ssize_t line;
    while ((line = getline(&buf, &n, puzzle_list)) >= 0) {
        if (line > 0 && buf[line - 1] == '\n')
            buf[line - 1] = '\0';
        struct puzzle_file *pf = parse_puzzle_file(buf);
        if (!pf) {
            fprintf(stderr, "couldn't parse puzzle at '%s'\n", buf);
            continue;
        }
        size_t last_slash = 0;
        for (size_t i = 0; buf[i]; ++i) {
            if (buf[i] == '/')
                last_slash = i + 1;
        }
        size_t last_dot = 0;
        for (size_t i = last_slash; buf[i]; ++i) {
            if (buf[i] == '.')
                last_dot = i;
        }
        if (last_dot > last_slash) {
            struct puzzle *puzzle = calloc(sizeof(struct puzzle), 1);
            puzzle->filename = malloc(last_dot - last_slash + 1);
            buf[last_dot] = '\0';
            memcpy(puzzle->filename, buf + last_slash, last_dot - last_slash + 1);
            puzzle->pf = pf;
            fprintf(stderr, "puzzle '%s' parsed\n", puzzle->filename);
            puzzle->next = puzzles;
            puzzles = puzzle;
        }
    }
    pclose(puzzle_list);

    int total_solutions = 0;
    int validated_solutions = 0;
    FILE *solution_list = popen("find test/solution -type f", "r");
    while ((line = getline(&buf, &n, solution_list)) >= 0) {
        if (line > 0 && buf[line - 1] == '\n')
            buf[line - 1] = '\0';
        struct solution_file *sf = parse_solution_file(buf);
        if (!sf) {
            fprintf(stderr, "couldn't parse solution at '%s'\n", buf);
            continue;
        }

        struct puzzle *puzzle = puzzles;
        while (puzzle && !byte_string_is(sf->puzzle, puzzle->filename))
            puzzle = puzzle->next;
        if (!puzzle) {
            fprintf(stderr, "couldn't find puzzle named '%.*s' for '%s'\n", (int)sf->puzzle.length, sf->puzzle.bytes, buf);
            free_solution_file(sf);
            continue;
        }

        struct solution solution = { 0 };
        struct board board = { 0 };
        if (!decode_solution(&solution, puzzle->pf, sf)) {
            free_solution_file(sf);
            continue;
        }
        total_solutions++;

        // set up the board.
        initial_setup(&solution, &board);

        // run the solution.
        while (board.cycle < 200000 && !board.complete) {
            cycle(&solution, &board);
            if (board.collision) {
                fprintf(stderr, "collision in '%s' at at %" PRId32 ", %" PRId32 ": %s\n", buf,
                 board.collision_location.u, board.collision_location.v,
                 board.collision_reason);
                break;
            }
        }
        if (sf->cycles != board.cycle) {
            fprintf(stderr, "cycle mismatch for '%s'\n", buf);
            fprintf(stderr, "solution file says cycle count is: %" PRIu32 "\n", sf->cycles);
            fprintf(stderr, "simulation says cycle count is: %" PRIu64 "\n", board.cycle);
        } else
            validated_solutions++;
        free_solution_file(sf);
    }
    pclose(puzzle_list);

    printf("%d / %d solutions validated!\n", validated_solutions, total_solutions);

    free(buf);
}
