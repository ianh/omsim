#ifndef OM_COLLISION_H
#define OM_COLLISION_H

#include "sim.h"

// this function also marks the area covered by any moving atoms.
bool collision(struct solution *solution, struct board *board, struct vector *collision_location);

#endif
