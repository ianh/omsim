#ifndef OM_SIM_H
#define OM_SIM_H

#include <stdbool.h>
#include <stdint.h>
#include <string.h>

typedef uint64_t atom;

#define NUMBER_OF_ATOM_TYPES 16
#define ATOM_OF_TYPE(type) (1ULL << ((type) + 1))

#define VALID (1ULL << 0)

// these shift amounts must match the atom bytes in the puzzle file format.
#define SALT (1ULL << 1)
#define AIR (1ULL << 2)
#define EARTH (1ULL << 3)
#define FIRE (1ULL << 4)
#define WATER (1ULL << 5)
#define QUICKSILVER (1ULL << 6)
#define GOLD (1ULL << 7)
#define SILVER (1ULL << 8)
#define COPPER (1ULL << 9)
#define IRON (1ULL << 10)
#define TIN (1ULL << 11)
#define LEAD (1ULL << 12)
#define VITAE (1ULL << 13)
#define MORS (1ULL << 14)
#define REPEATING_OUTPUT_PLACEHOLDER (1ULL << 15)
#define QUINTESSENCE (1ULL << 16)

// an output which accepts any molecule of the proper shape.  used for
// computation puzzles.
#define VARIABLE_OUTPUT (1ULL << 17)

// when removing this atom, replace it with one of the atoms in
// board->overlapped_atoms.  we can reuse the bit for VARIABLE_OUTPUT because
// that's only set on output atoms, whereas this bit never is.
#define OVERLAPS_ATOMS (1ULL << 17)

// conversion glyphs like animismus put this flag on their outputs until the
// second half-cycle.  it stops their outputs from being seen by other glyphs.
#define BEING_PRODUCED (1ULL << 18)

// conduits only transport atoms that have just been dropped.
#define BEING_DROPPED (1ULL << 19)

// is this atom part of a van berlo's wheel?
#define VAN_BERLO_ATOM (1ULL << 20)

// the motion of this atom is being tracked in a side table.  the flag and the
// atom's motion data will be cleared if the atom enters the board's bounding
// box.
#define IS_CHAIN_ATOM (1ULL << 21)

// is this atom being grabbed?  prevents output and consumption by glyphs.  the
// full 5-bit value is the number of times the atom has been grabbed (this is
// necessary to keep track of multiple simultaneous grabs).
#define GRABBED_ONCE (1ULL << 23)
#define GRABBED (0x1FULL * GRABBED_ONCE)

#define REMOVED (1ULL << 27)
// removed atoms can be marked with their movement index to detect movement in
// different directions.  these flags are only valid for removed atoms.
#define MOVED (1ULL << 28)
#define MOVEMENT_INDEX(movement_index) ((atom)(movement_index) << 29)
#define GET_MOVEMENT_INDEX(atom) ((atom) >> 29)
#define MAX_MOVEMENTS ((-1ULL) >> 29)

// offsets for the bits that indicate bonds.
#define RECENT_BOND 29 // also used to mark neighboring atoms in infinite products.
#define NORMAL_BOND 35
#define TRIPLEX_BOND_R 41
#define TRIPLEX_BOND_Y 47
#define TRIPLEX_BOND_K 53

// rotating can touch the 5 bits after the bonds (59-63), so make sure the
// following flags are clear before rotating a molecule.

// marks atoms inside a region of interest.  used for conduits.
#define CONDUIT_SHAPE (1ULL << 59)

// used for molecule flood fills.
#define VISITED (1ULL << 60)

#define BOND_LOW_BITS ((1ULL << RECENT_BOND) | (1ULL << NORMAL_BOND) | \
 (1ULL << TRIPLEX_BOND_R) | (1ULL << TRIPLEX_BOND_Y) | (1ULL << TRIPLEX_BOND_K))

#define ANY_ELEMENTAL (WATER | FIRE | EARTH | AIR)
#define ANY_METAL (LEAD | TIN | IRON | COPPER | SILVER | GOLD)
#define ANY_ATOM (SALT | ANY_ELEMENTAL | QUICKSILVER | ANY_METAL | \
 VITAE | MORS | QUINTESSENCE)
#define RECENT_BONDS (0x3FULL << RECENT_BOND)
#define NORMAL_BONDS (0x3FULL << NORMAL_BOND)
#define TRIPLEX_R_BONDS (0x3FULL << TRIPLEX_BOND_R)
#define TRIPLEX_Y_BONDS (0x3FULL << TRIPLEX_BOND_Y)
#define TRIPLEX_K_BONDS (0x3FULL << TRIPLEX_BOND_K)
#define TRIPLEX_BONDS (TRIPLEX_R_BONDS | TRIPLEX_Y_BONDS | TRIPLEX_K_BONDS)
#define ALL_BONDS (0x3FULL * BOND_LOW_BITS)

// TEMPORARY_FLAGS are reset in reset_temporary_flags().  the
// schedule_flag_reset_if_needed() function must be called whenever one of these
// flags is set.
#define TEMPORARY_FLAGS (BEING_DROPPED | RECENT_BONDS)

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

    CONDUIT = 1 << 13,

    ARM = 1 << 14,
    TWO_ARM = 1 << 15,
    THREE_ARM = 1 << 16,
    SIX_ARM = 1 << 17,
    PISTON = 1 << 18,

    VAN_BERLO = 1 << 20,

    // what is this arm grabbing?  each of the 6 possible grabbing directions
    // are tracked using a separate bit.
    GRABBING_LOW_BIT = 1 << 21,

    // is the arm itself grabbing?
    GRABBING = 1 << 27,

    // were the grabbed atoms moved this cycle?
    MOVED_GRABBED_ATOMS = 1 << 28,
};

static const enum mechanism_type GRABBING_EVERYTHING = 0x3FULL * GRABBING_LOW_BIT;

static const enum mechanism_type ANY_GLYPH = CALCIFICATION | ANIMISMUS |
 PROJECTION | DISPERSION | PURIFICATION | DUPLICATION | UNIFICATION | BONDING |
 UNBONDING | TRIPLEX_BONDING | MULTI_BONDING | DISPOSAL | EQUILIBRIUM | CONDUIT;
static const enum mechanism_type CONVERSION_GLYPH = ANIMISMUS | DISPERSION |
 PURIFICATION | UNIFICATION;

static const enum mechanism_type ANY_ARM = ARM | TWO_ARM | THREE_ARM | SIX_ARM |
 PISTON | VAN_BERLO;

struct vector {
    int32_t u;
    int32_t v;
};
static const struct vector zero_vector;
static inline bool vectors_equal(struct vector a, struct vector b)
{
    return a.u == b.u && a.v == b.v;
}
static inline uint32_t vector_hash(struct vector p)
{
    uint64_t hash = 0x542ddeaec5c75e0full;
    hash += 0xc0c594dd042705fbull * (uint32_t)p.u;
    hash += 0x93b45776c46130c1ull * (uint32_t)p.v;
    return hash >> 32;
}

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

    // the index of the conduit this glyph is associated with.  only meaningful
    // for mechanisms of type CONDUIT.
    uint32_t conduit_index;

    // used for arms that move on tracks.
    struct vector movement;

    // used for arms to track rotation for collision and "overclock" detection.
    int32_t arm_rotation;

    // used for visual loop detection.  either 0 or 1.
    int pivot_parity;
};

struct conduit {
    // the index of the glyph associated with this conduit.
    uint32_t glyph_index;

    // the index of the glyph associated with the other side of the conduit.
    uint32_t other_side_glyph_index;

    // the conduit id from the solution file.
    uint32_t id;

    // where the conduit's holes are relative to its glyph location.
    struct vector *positions;
    uint32_t number_of_positions;

    // the atoms in the conduit.  updated once at the beginning of each cycle
    // and again when the conduit glyph consumes.  when the conduit glyph
    // produces, these atoms appear at the other side of the conduit.  each
    // molecule in the conduit appears in order.
    // note that atom bonds are left unmodifed in the board's coordinate space,
    // while positions are transformed into the conduit glyph's coordinate
    // space.
    struct atom_at_position *atoms;

    // the number of atoms in each molecule in the conduit.
    uint32_t *molecule_lengths;
    uint32_t number_of_molecules;
};

enum input_output_type {
    INPUT = 1 << 0,
    SINGLE_OUTPUT = 1 << 1,
    REPEATING_OUTPUT = 1 << 2,

    // either single or repeating.
    OUTPUT = 3 << 1,

    // flag for inputs.
    BLOCKED = 1 << 3,

    // stop running and return INPUT_OUTPUT before the input spawns or the
    // output consumes.  you can use this flag to implement dynamic inputs and
    // outputs.
    INTERRUPT = 1 << 4,
};

// the number of times repeating outputs repeat.
#define REPEATING_OUTPUT_REPETITIONS 6

struct input_output {
    enum input_output_type type;

    struct atom_at_position *atoms;
    uint32_t number_of_atoms;

    // the original index of this input or output in the puzzle file.
    uint32_t puzzle_index;
    // the original index of this input or output in the solution file.
    uint32_t solution_index;

    // the number of times this output has consumed something.
    uint64_t number_of_outputs;

    // normally, a repeating output only goes up to six repetitions.  but we
    // keep counting after that so throughput can be measured in the limit.
    // that means extending the footprint of the output every time we hit a new
    // multiple of six.  this variable tracks the number of repetitions the
    // footprint currently contains.
    uint32_t number_of_repetitions;
    // how much should number_of_outputs go up for every matched repetition?
    uint32_t outputs_per_repetition;

    // in EXTEND_CHAIN mode, this field records the fastest rate at which a
    // matching output shape extends (in hexes per steady-state period).
    int32_t maximum_feed_rate;

    // these are the original atoms, before repetition is applied.  the
    // placeholder atom is guaranteed to be at the last position in the array.
    struct atom_at_position *original_atoms;
    uint32_t number_of_original_atoms;

    // where the repeated atoms will be attached to the repetition placeholder.
    struct vector repetition_origin;

    // bounding box information for repeating outputs.
    int32_t min_v;
    int32_t max_v;
    // these arrays each have length (max_v - min_v + 1).  row_min_u[i] is the
    // minimum output u coordinate for the row where v = min_v + i, and
    // row_max_u[i] is the corresponding maximum coordinate for that row.
    int32_t *row_min_u;
    int32_t *row_max_u;
};

struct solution {
    struct mechanism *glyphs;
    size_t number_of_glyphs;

    struct mechanism *arms;
    // array of instruction tape arrays, one per arm.
    char **arm_tape;
    // the lengths of these tape arrays.
    size_t *arm_tape_length;
    // the cycles on which to begin reading instructions from each tape.
    int64_t *arm_tape_start_cycle;
    size_t number_of_arms;

    // the maximum absolute value of all arm rotation amounts over all cycles.
    // used to detect "overclocking".
    uint32_t maximum_absolute_arm_rotation;

    // how many cycles until each tape loops back to the beginning.
    uint64_t tape_period;

    struct vector *track_positions;
    struct vector *track_plus_motions;
    struct vector *track_minus_motions;
    // track_table_size is a power of two, and the track_*** arrays form a hash
    // table keyed on track_positions.  a position of (INT32_MIN, INT32_MIN)
    // indicates an empty slot in the hash table.
    uint32_t track_table_size;

    // overlapping track is impossible to detect after a solution is decoded, so
    // the decoder sets this flag if it sees any.
    uint64_t track_self_overlap;

    struct conduit *conduits;
    size_t number_of_conduits;

    struct vector *cabinet_walls;
    size_t number_of_cabinet_walls;

    // whether it's an input or an output is determined by the type.
    struct input_output *inputs_and_outputs;
    size_t number_of_inputs_and_outputs;

    uint64_t target_number_of_outputs;

    // only atoms inside this bounding box are visible to the mechanisms of the
    // solution.
    int32_t min_visible_u;
    int32_t max_visible_u;
    int32_t min_visible_v;
    int32_t max_visible_v;
};
enum movement_type {
    SWING_MOVEMENT = 0,
    PIVOT_MOVEMENT = 1,
    PISTON_MOVEMENT = 2,
    TRACK_MOVEMENT = 3,

    IS_PISTON = 4,
};
struct movement {
    // the number of atoms involved in this movement.
    size_t number_of_atoms;
    size_t first_atom_index;

    uint32_t type;

    struct vector base;
    struct vector grabber_offset; // without rotations applied.
    struct vector absolute_grab_position;
    struct vector translation;

    int32_t base_rotation;
    int32_t rotation;

    int piston_extension;

    // linked list of chain_atom structs corresponding to atoms
    // involved in this movement.
    uint32_t first_chain_atom;
};

struct movement_list {
    struct movement *movements;
    size_t capacity;
    size_t length;
    size_t cursor;
};
struct moving_atoms {
    struct atom_at_position *atoms_at_positions;
    size_t capacity;
    size_t length;
    size_t cursor;
};
struct marked_positions {
    struct vector *positions;
    size_t capacity;
    size_t length;
};

// chain atoms participate in waste chain / polymer throughput detection.
#define CHAIN_ATOM_ROTATION 7u
#define CHAIN_ATOM_IN_REPEATING_SEGMENT (1u << 3)
#define CHAIN_ATOM_SWING_SEXTANTS_SHIFT 4
#define CHAIN_ATOM_SWING_SEXTANTS (63u << CHAIN_ATOM_SWING_SEXTANTS_SHIFT)
struct chain_atom {
    uint32_t *prev_in_list;
    uint32_t next_in_list;
    uint32_t flags;
    struct vector current_position;
    struct vector original_position;
    uint32_t area_direction;
};
// in the DISCOVER_CHAIN mode, the original position is left fixed so that the
// overall motion can be discovered.  in the EXTEND_CHAIN mode, the original
// position is carried along with movements to track the relative translation
// throughout the steady state period.
enum chain_mode {
    DISCOVER_CHAIN,
    EXTEND_CHAIN,
};
// used for chain atom area.
enum growth_order {
    GROWTH_NONE,
    GROWTH_LINEAR,
    GROWTH_QUADRATIC,
};
#define GRID_ARRAY_MIN (-32)
#define GRID_ARRAY_MAX 32
#define GRID_ARRAY_PREFIX ((GRID_ARRAY_MAX-GRID_ARRAY_MIN)*(GRID_ARRAY_MAX-GRID_ARRAY_MIN))
#define GRID_CAPACITY(grid) (GRID_ARRAY_PREFIX + (grid).hash_capacity)
#define BOARD_CAPACITY(board) (GRID_CAPACITY((board)->grid))
struct atom_grid {
    struct atom_at_position *atoms_at_positions;
    uint32_t hash_capacity;
    struct board *board;
};
struct linear_area_direction {
    struct vector direction;
    struct atom_grid footprint_at_infinity;
};
struct board {
    struct atom_grid grid;

    uint32_t area;

    // atoms that need their flags reset.
    atom **flag_reset;
    uint32_t flag_reset_capacity;
    uint32_t flag_reset_length;

    uint64_t cycle;
    int half_cycle;

    // set whenever run() returns INPUT_OUTPUT.
    size_t active_input_or_output;

    struct movement_list movements;
    struct moving_atoms moving_atoms;

    struct atom_at_position *overlapped_atoms;
    uint32_t number_of_overlapped_atoms;
    uint32_t overlapped_atoms_capacity;
    const char *overlapped_reason;

    // used for checking infinite products.
    struct marked_positions marked;

    bool collision;
    struct vector collision_location;
    const char *collision_reason;

    uint32_t atom_grabs[NUMBER_OF_ATOM_TYPES];

    // a bitmask of which output indexes cause the solution to fail if a wrong
    // output of the correct shape is dropped onto them.
    uint64_t fails_on_wrong_output_mask;
    // a bitmask of which output indexes cause the solution to fail if a an
    // output of the correct footprint but with wrong bonds is dropped onto them.
    uint64_t fails_on_wrong_output_bonds_mask;
    // the index of the wrong output in the solution file (or SIZE_MAX).
    size_t wrong_output_index;

    // how many hexes overlap one another (other than arms and track)?
    uint64_t overlap;

    enum chain_mode chain_mode;
    // this is the side table for atoms with IS_CHAIN_ATOM set.
    struct chain_atom *chain_atoms;
    uint32_t number_of_chain_atoms;
    // a hash table mapping atom positions to their index in chain_atoms.
    uint32_t *chain_atom_table;
    uint32_t chain_atom_table_size;
    // set if a chain atom re-enters the box while chain_mode is EXTEND_CHAIN.
    bool chain_will_become_visible;

    enum growth_order area_growth_order;
    struct linear_area_direction *area_directions;
    uint32_t number_of_area_directions;

    // records each cycle the output count increases toward completion.
    // a cycle on which the output count increased multiple times will appear
    // that number of times in the list.
    uint64_t *output_cycles;
    uint64_t number_of_output_cycles;
    uint64_t output_cycles_capacity;

    // how many times a collision check has been performed.
    uint64_t collision_checks;
    // the limit after which execution will stop.
    uint64_t collision_check_limit;

    // did the solution complete yet?
    bool complete;
};

// initial_board_size is a sizing hint to the hash table; zero means 'no hint'.
void initial_setup(struct solution *solution, struct board *board,
 uint32_t intial_board_size);

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

void insert_atom(struct board *board, struct vector query, atom atom, const char *collision_reason);
atom *lookup_atom(struct board *board, struct vector query);
atom *lookup_atom_in_grid(struct atom_grid *grid, struct vector query);
struct atom_at_position *lookup_atom_at_position(struct atom_grid *grid, struct vector query);

// returns the atom value at the point (for collision detection).
atom mark_used_area(struct board *board, struct vector point);

bool lookup_track(struct solution *solution, struct vector query, uint32_t *index);

uint32_t used_area(struct board *board);

static inline bool position_may_be_visible_to_solution(struct solution *solution, struct vector p)
{
    return p.u >= solution->min_visible_u && p.u <= solution->max_visible_u &&
        p.v >= solution->min_visible_v && p.v <= solution->max_visible_v;
}
void add_chain_atom_to_table(struct board *board, uint32_t chain_atom_index);
uint32_t lookup_chain_atom(struct board *board, struct vector query);
void move_chain_atom_to_list(struct board *board, uint32_t chain_atom_index, uint32_t *list);

// used during decoding.
bool repeat_molecule(struct input_output *io, uint32_t number_of_repetitions,
 const char **error);

// the origin is always the last vector in the footprint of a glyph.
const struct vector *glyph_footprint(enum mechanism_type type);

// geometric helper functions.
struct vector u_offset_for_direction(int direction);
struct vector v_offset_for_direction(int direction);
int direction_for_offset(struct vector d);
int angular_distance_between_grabbers(enum mechanism_type);
struct vector mechanism_relative_position(struct mechanism m, int32_t du, int32_t dv, int32_t w);
atom bond_direction(struct mechanism m, int32_t du, int32_t dv);

#endif
