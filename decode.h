#ifndef OM_DECODE_H
#define OM_DECODE_H

#include <stdbool.h>
#include <stdint.h>

#define GLYPH_OUTSIDE_CABINET (1UL << 0)
#define ARM_OUTSIDE_CABINET (1UL << 1)
#define ARM_REACHES_ACROSS_WALL (1UL << 2)
#define INPUT_OUTPUT_OUTSIDE_CABINET (1UL << 3)
#define TRACK_OUTSIDE_CABINET (1UL << 4)
#define ISOLATION_VIOLATED (1UL << 5)
#define CONDUIT_IN_WRONG_CABINET (1UL << 6)
#define CONDUIT_ALTERED (1UL << 7)
#define CONDUIT_IN_FREESPACE (1UL << 8)

struct puzzle_file;
struct solution;
struct solution_file;

bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf, const char **error);

uint64_t solution_file_cost(struct solution_file *sf);
uint64_t solution_instructions(struct solution *solution);

#endif
