#ifndef OM_DECODE_H
#define OM_DECODE_H

#include <stdbool.h>
#include <stdint.h>

#define CONDUIT_VIOLATED (1UL << 0)
#define GLYPH_VIOLATED (1UL << 1)
#define ARM_VIOLATED (1UL << 2)
#define CROSS_WALL_VIOLATED (1UL << 3)
#define INPUT_OUTPUT_VIOLATED (1UL << 4)
#define TRACK_VIOLATED (1UL << 5)
#define ISOLATION_VIOLATED (1UL << 6)

struct puzzle_file;
struct solution;
struct solution_file;

bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf, const char **error);

uint64_t solution_file_cost(struct solution_file *sf);
uint64_t solution_instructions(struct solution *solution);

#endif
