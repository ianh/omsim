#ifndef OM_DECODE_H
#define OM_DECODE_H

#include <stdbool.h>

struct puzzle_file;
struct solution;
struct solution_file;

bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf);

#endif
