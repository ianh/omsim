#ifndef OM_DECODE_H
#define OM_DECODE_H

#include <stdbool.h>
#include <stdint.h>

struct puzzle_file;
struct solution;
struct solution_file;

bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf, const char **error);

uint64_t solution_file_cost(struct solution_file *sf);
uint64_t solution_instructions(struct solution *solution);

#endif
