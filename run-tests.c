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

#define AREA_TOLERANCE 0.008

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
    FILE *solution_list = popen("find test/solution -type f -name *.solution", "r");
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
        const char *error;
        if (!decode_solution(&solution, puzzle->pf, sf, &error)) {
            fprintf(stderr, "error in '%s': %s\n", buf, error);
            free_solution_file(sf);
            continue;
        }
        total_solutions++;

        uint64_t cost = solution_file_cost(sf);
        if (sf->cost != cost) {
            fprintf(stderr, "cost mismatch for '%s'\n", buf);
            fprintf(stderr, "solution file says cost is: %" PRIu32 "\n", sf->cost);
            fprintf(stderr, "adding up its parts, the cost is: %" PRIu64 "\n", cost);
            goto fail;
        }

        uint64_t instructions = solution_instructions(&solution);
        if (sf->instructions != instructions) {
            fprintf(stderr, "instructions mismatch for '%s'\n", buf);
            fprintf(stderr, "solution file says instruction count is: %" PRIu32 "\n", sf->instructions);
            fprintf(stderr, "counting instructions says instruction count is: %" PRIu64 "\n", instructions);
            goto fail;
        }

        // set up the board.
        initial_setup(&solution, &board, sf->area);

        // run the solution.
        while (board.cycle < 200000 && !board.complete) {
            cycle(&solution, &board);
            if (board.collision) {
                fprintf(stderr, "collision in '%s' at %" PRId32 ", %" PRId32 ": %s\n", buf,
                 board.collision_location.u, board.collision_location.v,
                 board.collision_reason);
                break;
            }
        }
        if (sf->cycles != board.cycle) {
            fprintf(stderr, "cycle mismatch for '%s'\n", buf);
            fprintf(stderr, "solution file says cycle count is: %" PRIu32 "\n", sf->cycles);
            fprintf(stderr, "simulation says cycle count is: %" PRIu64 "\n", board.cycle);
            goto fail;
        }
        uint32_t area = used_area(&board);
        if (!puzzle->pf->production_info && sf->area != area) {
            fprintf(stderr, "area mismatch for '%s'\n", buf);
            fprintf(stderr, "solution file says area is: %" PRIu32 "\n", sf->area);
            fprintf(stderr, "simulation says area is: %" PRIu32 "\n", area);
            goto fail;
        }
        validated_solutions++;
    fail:
        destroy(&solution, &board);
        free_solution_file(sf);
    }
    pclose(solution_list);

    fprintf(stderr, "%d / %d solutions validated!\n", validated_solutions, total_solutions);

    free(buf);

    struct puzzle *puzzle = puzzles;
    while (puzzle) {
        struct puzzle *next = puzzle->next;
        free_puzzle_file(puzzle->pf);
        free(puzzle->filename);
        free(puzzle);
        puzzle = next;
    }
}
