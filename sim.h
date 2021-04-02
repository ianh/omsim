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

// this flag is necessary to simulate overlapping glyphs correctly.
static const atom UNBONDED = 1ULL << 17;

// conversion glyphs like animismus put this flag on their outputs until the
// second half-cycle.  it stops their outputs from being seen by other glyphs.
static const atom BEING_PRODUCED = 1ULL << 18;

// is this atom part of a van berlo's wheel?
static const atom VAN_BERLO_ATOM = 1ULL << 19;

// is this atom being grabbed?  prevents output and consumption by glyphs.  the
// full 5-bit value is the number of times the atom has been grabbed (this is
// necessary to keep track of multiple simultaneous grabs).
static const atom GRABBED_ONCE = 1ULL << 20;
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

struct mechanism {
    enum mechanism_type type;

    struct vector position;

    // direction (two basis vectors).  includes length for arms.
    struct vector direction_u;
    struct vector direction_v;
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

    // input_atoms contains the atoms from the first input molecule, then the
    // atoms from the second input molecule, and so on.
    atom *input_atoms;
    struct vector *input_atom_positions;
    uint32_t *input_molecule_lengths;
    bool *input_molecule_is_blocked;
    // the original index of this molecule in the puzzle file.
    uint32_t *input_molecule_puzzle_index;
    uint32_t number_of_input_molecules;

    atom *output_atoms;
    struct vector *output_atom_positions;
    uint32_t *output_molecule_lengths;
    uint64_t *output_molecule_number_of_outputs;
    uint32_t number_of_output_molecules;

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
    struct {
        struct vector position;
        atom atom;
    } *atoms_at_positions;

    uint32_t capacity;
    uint32_t used;
    uint32_t removed;

    uint64_t cycle;
    uint32_t half_cycle;

    struct movement_list movements;

    bool collision;
    struct vector collision_location;
    const char *collision_reason;

    bool complete;
};

struct puzzle_file;
struct solution_file;

bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf);
void initial_setup(struct solution *solution, struct board *board);
void cycle(struct solution *solution, struct board *board);

#endif
