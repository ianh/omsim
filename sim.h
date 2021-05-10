#ifndef OM_SIM_H
#define OM_SIM_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

typedef uint64_t atom;

// these shift amounts must match the atom bytes in the puzzle file format.
static const atom SALT = 1ULL << 1;
static const atom AIR = 1ULL << 2;
static const atom EARTH = 1ULL << 3;
static const atom FIRE = 1ULL << 4;
static const atom WATER = 1ULL << 5;
static const atom QUICKSILVER = 1ULL << 6;
static const atom GOLD = 1ULL << 7;
static const atom SILVER = 1ULL << 8;
static const atom COPPER = 1ULL << 9;
static const atom IRON = 1ULL << 10;
static const atom TIN = 1ULL << 11;
static const atom LEAD = 1ULL << 12;
static const atom VITAE = 1ULL << 13;
static const atom MORS = 1ULL << 14;
// this is commented out because it's unsupported (for now).
// static const atom REPEATING_OUTPUT = 1ULL << 15;
static const atom QUINTESSENCE = 1ULL << 16;

// an output which accepts any molecule of the proper shape.  used for
// computation puzzles.
static const atom VARIABLE_OUTPUT = 1ULL << 17;

// this flag is necessary to simulate overlapping glyphs correctly.
static const atom UNBONDED = 1ULL << 18;

// conversion glyphs like animismus put this flag on their outputs until the
// second half-cycle.  it stops their outputs from being seen by other glyphs.
static const atom BEING_PRODUCED = 1ULL << 19;

// is this atom part of a van berlo's wheel?
static const atom VAN_BERLO_ATOM = 1ULL << 20;

// is this atom being grabbed?  prevents output and consumption by glyphs.  the
// full 5-bit value is the number of times the atom has been grabbed (this is
// necessary to keep track of multiple simultaneous grabs).
static const atom GRABBED_ONCE = 1ULL << 25;
static const atom GRABBED = 0x1FULL * GRABBED_ONCE;

static const atom VALID = 1ULL << 30;
static const atom REMOVED = 1ULL << 31;

// offsets for the bits that indicate bonds.
static const int NORMAL_BOND = 32;
static const int TRIPLEX_BOND_R = 38;
static const int TRIPLEX_BOND_Y = 44;
static const int TRIPLEX_BOND_K = 50;
// ensure that the 5 bits after aren't used for anything else, since rotating
// atoms can touch these bits temporarily.

static const atom BOND_LOW_BITS = (1ULL << NORMAL_BOND) |
 (1ULL << TRIPLEX_BOND_R) | (1ULL << TRIPLEX_BOND_Y) |
 (1ULL << TRIPLEX_BOND_K);

static const atom ANY_ELEMENTAL = WATER | FIRE | EARTH | AIR;
static const atom ANY_METAL = LEAD | TIN | IRON | COPPER | SILVER | GOLD;
static const atom ANY_ATOM = SALT | ANY_ELEMENTAL | QUICKSILVER | ANY_METAL |
 VITAE | MORS | QUINTESSENCE;
static const atom NORMAL_BONDS = 0x3FULL << NORMAL_BOND;
static const atom TRIPLEX_R_BONDS = 0x3FULL << TRIPLEX_BOND_R;
static const atom TRIPLEX_Y_BONDS = 0x3FULL << TRIPLEX_BOND_Y;
static const atom TRIPLEX_K_BONDS = 0x3FULL << TRIPLEX_BOND_K;
static const atom TRIPLEX_BONDS = TRIPLEX_R_BONDS | TRIPLEX_Y_BONDS | TRIPLEX_K_BONDS;
static const atom ALL_BONDS = NORMAL_BONDS | TRIPLEX_BONDS;

enum mechanism_type {
    NO_MECHANISM,

    CALCIFICATION = 1 << 0,
    ANIMISMUS = 1 << 1,
    PROJECTION = 1 << 2,
    DISPERSION = 1 << 3,
    PURIFICATION = 1 << 4,
    DUPLICATION = 1 << 5,
    UNIFICATION = 1 << 6,
    BONDING = 1 << 7,
    UNBONDING = 1 << 8,
    TRIPLEX_BONDING = 1 << 9,
    MULTI_BONDING = 1 << 10,
    DISPOSAL = 1 << 11,
    EQUILIBRIUM = 1 << 12,

    ARM = 1 << 13,
    TWO_ARM = 1 << 14,
    THREE_ARM = 1 << 15,
    SIX_ARM = 1 << 16,
    PISTON = 1 << 17,

    VAN_BERLO = 1 << 20,

    // what is this arm grabbing?  each of the 6 possible grabbing directions
    // are tracked using a separate bit.
    GRABBING_LOW_BIT = 1 << 21,

    // is the arm itself grabbing?
    GRABBING = 1 << 27,
};

static const enum mechanism_type GRABBING_EVERYTHING = 0x3FULL * GRABBING_LOW_BIT;

static const enum mechanism_type ANY_GLYPH = CALCIFICATION | ANIMISMUS |
 PROJECTION | DISPERSION | PURIFICATION | DUPLICATION | UNIFICATION | BONDING |
 UNBONDING | TRIPLEX_BONDING | MULTI_BONDING | DISPOSAL | EQUILIBRIUM;
static const enum mechanism_type CONVERSION_GLYPH = ANIMISMUS | DISPERSION |
 PURIFICATION | UNIFICATION;

static const enum mechanism_type ANY_ARM = ARM | TWO_ARM | THREE_ARM | SIX_ARM |
 PISTON | VAN_BERLO;

struct vector {
    int32_t u;
    int32_t v;
};

struct atom_at_position {
    struct vector position;
    atom atom;
};

struct mechanism {
    enum mechanism_type type;

    struct vector position;

    // direction (two basis vectors).  includes length for arms.
    struct vector direction_u;
    struct vector direction_v;
};

enum input_output_type {
    INPUT = 1 << 0,
    OUTPUT = 1 << 1,

    // flag for inputs.
    BLOCKED = 1 << 2,

    // stop running and return INPUT_OUTPUT before the input spawns or the
    // output consumes.  you can use this flag to implement dynamic inputs and
    // outputs.
    INTERRUPT = 1 << 3,
};

struct input_output {
    enum input_output_type type;

    struct atom_at_position *atoms;
    uint32_t number_of_atoms;

    // the original index of this input or output in the puzzle file.
    uint32_t puzzle_index;

    // output-specific fields.
    uint64_t number_of_outputs;
};

struct solution {
    struct mechanism *glyphs;
    size_t number_of_glyphs;

    struct mechanism *arms;
    // array of arrays.
    char **arm_tape;
    size_t *arm_tape_length;
    uint64_t *arm_tape_start_cycle;
    size_t number_of_arms;

    // how many cycles until each tape loops back to the beginning.
    uint64_t tape_period;

    struct vector *track_positions;
    struct vector *track_plus_motions;
    struct vector *track_minus_motions;
    // track_table_size is a power of two, and the track_*** arrays form a hash
    // table keyed on track_positions.  a position of (INT32_MIN, INT32_MIN)
    // indicates an empty slot in the hash table.
    uint32_t track_table_size;

    // whether it's an input or an output is determined by the type.
    struct input_output *inputs_and_outputs;
    size_t number_of_inputs_and_outputs;

    uint64_t target_number_of_outputs;
};

struct movement {
    // the starting position of the atom.
    struct vector position;

    // the atom being moved.
    atom atom;

    // for rotations, the base is the point around which the atom rotates.
    // for translations, the base is the direction of translation.
    struct vector base;
    int rotation;
    int translation;
};
struct movement_list {
    struct movement *elements;
    size_t capacity;
    size_t length;
    size_t cursor;
};
struct board {
    struct atom_at_position *atoms_at_positions;

    uint32_t capacity;
    uint32_t used;
    uint32_t removed;

    uint64_t cycle;
    int half_cycle;

    // set whenever run() returns INPUT_OUTPUT.
    size_t active_input_or_output;

    struct movement_list movements;

    bool collision;
    struct vector collision_location;
    const char *collision_reason;

    bool complete;
};

void initial_setup(struct solution *solution, struct board *board);

enum run_result {
    FINISHED_CYCLE,

    // the active input/output is available in board->active_input_or_output.
    INPUT_OUTPUT,
};
enum run_result run(struct solution *solution, struct board *board);

// shortcut function if you aren't using dynamic inputs or outputs.
static inline void cycle(struct solution *solution, struct board *board)
{
    while (run(solution, board) != FINISHED_CYCLE);
}

// free memory associated with the solution and board.
void destroy(struct solution *solution, struct board *board);

atom *lookup_atom(struct board *board, struct vector query);
bool lookup_track(struct solution *solution, struct vector query, uint32_t *index);

// geometric helper functions.
struct vector u_offset_for_direction(int direction);
struct vector v_offset_for_direction(int direction);
int direction_for_offset(struct vector d);
struct vector mechanism_relative_position(struct mechanism m, int32_t du, int32_t dv, int32_t w);
atom bond_direction(struct mechanism m, int32_t du, int32_t dv);

#endif
