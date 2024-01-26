#ifndef OM_STEADY_STATE_H
#define OM_STEADY_STATE_H

#include "sim.h"

enum eventual_behavior {
    EVENTUALLY_ENTERS_STEADY_STATE,
    EVENTUALLY_STOPS_RUNNING,
    EVENTUALLY_REACHES_CYCLE_LIMIT,
};

enum growth_order {
    GROWTH_NONE,
    GROWTH_LINEAR,
    GROWTH_QUADRATIC,
};

struct steady_state {
    // what type of overall behavior does the solution have in the limit?
    enum eventual_behavior eventual_behavior;

    // how many cycles does it take to loop in the steady state?
    uint64_t number_of_cycles;

    // how many outputs are produced during this loop?
    uint64_t number_of_outputs;

    // after which cycle do outputs start to repeat?
    uint64_t outputs_repeat_after_cycle;

    // is there a visual difference between the start and end of the loop due to
    // arm grabbers pivoting?
    bool pivot_parity;
};

struct steady_state run_until_steady_state(struct solution *solution, struct board *board, uint64_t cycle_limit);

#endif
