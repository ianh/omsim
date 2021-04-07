#include "parse.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

// hash table functions -- see appendix.
static bool lookup_track(struct solution *solution, struct vector query, uint32_t *index);
static atom *lookup_atom(struct board *board, struct vector query);
static atom *insert_atom(struct board *board, struct vector query, const char *collision_reason);
static void ensure_capacity(struct board *board, uint32_t potential_insertions);

static struct vector mechanism_relative_position(struct mechanism m, int32_t du, int32_t dv, int32_t w)
{
    return (struct vector){
        m.direction_u.u * du + m.direction_v.u * dv + m.position.u * w,
        m.direction_u.v * du + m.direction_v.v * dv + m.position.v * w,
    };
}

static void report_collision(struct board *board, struct vector p, const char *reason)
{
    // only report the first collision (after that, all bets are off).
    if (board->collision)
        return;
    board->collision = true;
    board->collision_location = p;
    board->collision_reason = reason;
}

static bool conversion_output(struct board *board, struct mechanism m, int32_t du, int32_t dv)
{
    assert(m.type & CONVERSION_GLYPH);
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    atom *output = lookup_atom(board, pos);
    if ((*output & VALID) && !(*output & REMOVED)) {
        if (board->half_cycle == 2) {
            // conversion glyph outputs appear in the second half-cycle.
            *output &= ~BEING_PRODUCED;
        } else if (*output & BEING_PRODUCED) {
            report_collision(board, pos, "two conversion glyphs outputting to the same point");
            return true;
        } else if (*output & VAN_BERLO_ATOM) {
            report_collision(board, pos, "conversion glyph output overlaps with van berlo's wheel");
            return true;
        }
        // xx add hitbox
        return false;
    }
    return true;
}

static void produce_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv, atom a)
{
    assert(m.type & CONVERSION_GLYPH);
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    *insert_atom(board, pos, "conversion glyph output") = VALID | BEING_PRODUCED | a;
}

static atom *get_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv)
{
    static const atom empty;
    // conversion glyphs don't consume any inputs in the second half-cycle.
    if ((m.type & CONVERSION_GLYPH) && board->half_cycle == 2)
        return (atom *)&empty;
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    atom *a = lookup_atom(board, pos);
    if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED))
        return (atom *)&empty;
    // only duplication glyphs can see the van berlo's wheel.
    if (!(m.type & (DUPLICATION | VAN_BERLO)) && (*a & VAN_BERLO_ATOM))
        return (atom *)&empty;
    // conversion glyphs can't see bonded or grabbed inputs.
    if ((m.type & CONVERSION_GLYPH) && ((*a & ALL_BONDS) || !(*a & UNBONDED) || (*a & GRABBED)))
        return (atom *)&empty;
    return a;
}

static inline void remove_atom(struct board *board, atom *a)
{
    assert(!(*a & REMOVED));
    board->used--;
    board->removed++;
    *a |= REMOVED;
}

static inline void transform_atom(atom *a, atom new_type)
{
    *a &= ~ANY_ATOM;
    *a |= new_type;
}

static int direction_for_offset(struct vector d)
{
    if (d.u == 0 && d.v == 1)
        return 0;
    else if (d.u == 1 && d.v == 0)
        return 1;
    else if (d.u == 1 && d.v == -1)
        return 2;
    else if (d.u == 0 && d.v == -1)
        return 3;
    else if (d.u == -1 && d.v == 0)
        return 4;
    else if (d.u == -1 && d.v == 1)
        return 5;
    return -1;
}

static atom bond_direction(struct mechanism m, int32_t du, int32_t dv)
{
    atom base = BOND_LOW_BITS;
    int dir = direction_for_offset(mechanism_relative_position(m, du, dv, 0));
    if (dir < 0) {
        fprintf(stderr, "internal error: non-orthonormal bond direction (scaled bonder?)\n");
        return 0;
    } else
        return base << dir;
}

// add the proper bond from the `bond` mask unless a bond from the `unless` mask
// already exists in that direction. 
static void add_bond(struct mechanism m, atom *a, atom *b, int32_t u, int32_t v, atom bond, atom unless)
{
    atom ab = bond_direction(m, u, v);
    atom ba = bond_direction(m, -u, -v);
    if ((*a & ab & unless) || (*b & ba & unless))
        return;
    *a |= ab & bond;
    *b |= ba & bond;
}

static int normalize_direction(int direction)
{
    direction %= 6;
    return direction < 0 ? direction + 6 : direction;
}

static struct vector v_offset_for_direction(int direction)
{
    switch (normalize_direction(direction)) {
    case 0: return (struct vector){0, 1};
    case 1: return (struct vector){1, 0};
    case 2: return (struct vector){1, -1};
    case 3: return (struct vector){0, -1};
    case 4: return (struct vector){-1, 0};
    case 5: return (struct vector){-1, 1};
    default: abort();
    }
}

static struct vector u_offset_for_direction(int direction)
{
    return v_offset_for_direction(direction + 1);
}

static void apply_glyphs(struct solution *solution, struct board *board)
{
    size_t n = solution->number_of_glyphs;
    for (size_t i = 0; i < n; ++i) {
        // at most 4 new atoms can be inserted by a single glyph (dispersion).
        // ensure there's enough space in the hash table for these new atoms.
        ensure_capacity(board, 4);
        struct mechanism m = solution->glyphs[i];
        switch (m.type & ANY_GLYPH) {
        case CALCIFICATION: {
            atom *a = get_atom(board, m, 0, 0);
            if (*a & ANY_ELEMENTAL)
                transform_atom(a, SALT);
            break;
        }
        case ANIMISMUS: {
            atom *a = get_atom(board, m, 0, 0);
            atom *b = get_atom(board, m, 0, 1);
            bool c = conversion_output(board, m, 1, 0);
            bool d = conversion_output(board, m, -1, 1);
            if (c && d && (*a & *b & SALT)) {
                remove_atom(board, a);
                remove_atom(board, b);
                produce_atom(board, m, 1, 0, VITAE);
                produce_atom(board, m, -1, 1, MORS);
            }
            break;
        }
        case PROJECTION: {
            atom *q = get_atom(board, m, 0, 0);
            atom *a = get_atom(board, m, 0, 1);
            atom metal = *a & ANY_METAL & ~GOLD;
            if (metal && (*q & UNBONDED) && !(*q & ALL_BONDS) && !(*q & GRABBED) && (*q & QUICKSILVER)) {
                remove_atom(board, q);
                transform_atom(a, metal >> 1);
            }
            break;
        }
        case DISPERSION: {
            atom *a = get_atom(board, m, 0, 0);
            bool b = conversion_output(board, m, 0, 1);
            bool c = conversion_output(board, m, -1, 1);
            bool d = conversion_output(board, m, -1, 0);
            bool e = conversion_output(board, m, 0, -1);
            if (b && c && d && e && (*a & QUINTESSENCE)) {
                remove_atom(board, a);
                produce_atom(board, m, 0, 1, EARTH);
                produce_atom(board, m, -1, 1, WATER);
                produce_atom(board, m, -1, 0, FIRE);
                produce_atom(board, m, 0, -1, AIR);
            }
            break;
        }
        case PURIFICATION: {
            atom *a = get_atom(board, m, 0, 0);
            atom *b = get_atom(board, m, 0, 1);
            bool c = conversion_output(board, m, 1, 0);
            atom metal = *a & *b & ANY_METAL & ~GOLD;
            if (c && metal) {
                remove_atom(board, a);
                remove_atom(board, b);
                produce_atom(board, m, 1, 0, metal >> 1);
            }
            break;
        }
        case DUPLICATION: {
            atom *a = get_atom(board, m, 0, 0);
            atom *b = get_atom(board, m, 0, 1);
            atom elemental = *a & ANY_ELEMENTAL;
            if (elemental && (*b & SALT) && !(*b & VAN_BERLO_ATOM))
                transform_atom(b, elemental);
            break;
        }
        case UNIFICATION: {
            atom *a = get_atom(board, m, 1, 0);
            atom *b = get_atom(board, m, 1, -1);
            atom *c = get_atom(board, m, -1, 0);
            atom *d = get_atom(board, m, -1, 1);
            bool e = conversion_output(board, m, 0, 0);
            if (e && ((*a | *b | *c | *d) & ANY_ELEMENTAL) == ANY_ELEMENTAL) {
                remove_atom(board, a);
                remove_atom(board, b);
                remove_atom(board, c);
                remove_atom(board, d);
                produce_atom(board, m, 0, 0, QUINTESSENCE);
            }
            break;
        }
        case BONDING: {
            atom *a = get_atom(board, m, 0, 0);
            atom *b = get_atom(board, m, 0, 1);
            if (*a && *b)
                add_bond(m, a, b, 0, 1, NORMAL_BONDS, TRIPLEX_BONDS);
            break;
        }
        case UNBONDING: {
            atom *a = get_atom(board, m, 0, 0);
            atom *b = get_atom(board, m, 0, 1);
            if (*a && *b) {
                *a &= ~bond_direction(m, 0, 1);
                *b &= ~bond_direction(m, 0, -1);
            }
            break;
        }
        case TRIPLEX_BONDING: {
            atom *ky = get_atom(board, m, 0, 0);
            atom *yr = get_atom(board, m, 1, 0);
            atom *rk = get_atom(board, m, 0, 1);
            if (*ky && *yr && (*ky & *yr & FIRE))
                add_bond(m, ky, yr, 1, 0, TRIPLEX_Y_BONDS, NORMAL_BONDS);
            if (*yr && *rk && (*yr & *rk & FIRE))
                add_bond(m, yr, rk, -1, 1, TRIPLEX_R_BONDS, NORMAL_BONDS);
            if (*rk && *ky && (*rk & *ky & FIRE))
                add_bond(m, rk, ky, 0, -1, TRIPLEX_K_BONDS, NORMAL_BONDS);
            break;
        }
        case MULTI_BONDING: {
            atom *center = get_atom(board, m, 0, 0);
            atom *a = get_atom(board, m, 0, 1);
            atom *b = get_atom(board, m, -1, 0);
            atom *c = get_atom(board, m, 1, -1);
            if (*center && *a)
                add_bond(m, a, center, 0, -1, NORMAL_BONDS, TRIPLEX_BONDS);
            if (*center && *b)
                add_bond(m, b, center, 1, 0, NORMAL_BONDS, TRIPLEX_BONDS);
            if (*center && *c)
                add_bond(m, c, center, -1, 1, NORMAL_BONDS, TRIPLEX_BONDS);
            break;
        }
        case DISPOSAL: {
            atom *a = get_atom(board, m, 0, 0);
            if (*a && (*a & UNBONDED) && !(*a & ALL_BONDS) && !(*a & GRABBED))
                remove_atom(board, a);
            break;
        }
        case EQUILIBRIUM:
        default:
            break;
        }
    }
}

static void adjust_axis_magnitude(struct vector *p, int32_t delta)
{
    if (p->u > 0)
        p->u += delta;
    if (p->u < 0)
        p->u -= delta;
    if (p->v > 0)
        p->v += delta;
    if (p->v < 0)
        p->v -= delta;
}

static struct vector normalize_axis(struct vector p)
{
    if (p.u > 0)
        p.u = 1;
    if (p.u < 0)
        p.u = -1;
    if (p.v > 0)
        p.v = 1;
    if (p.v < 0)
        p.v = -1;
    return p;
}

static void enqueue_movement(struct board *board, struct movement m)
{
    size_t capacity = board->movements.capacity;
    while (board->movements.length >= capacity)
        capacity = 4 * (capacity + 12) / 3;
    if (capacity != board->movements.capacity) {
        struct movement *elems = realloc(board->movements.elements,
         sizeof(struct movement) * capacity);
        if (!elems)
            abort();
        board->movements.elements = elems;
        board->movements.capacity = capacity;
    }
    board->movements.elements[board->movements.length++] = m;
}

static void move_atoms(struct board *board, atom *a, struct vector position, struct vector base, int rotation, int translation)
{
    enqueue_movement(board, (struct movement){
        .position = position,
        .atom = *a,
        .base = base,
        .rotation = rotation,
        .translation = translation,
    });
    remove_atom(board, a);
    // do a breadth-first search over the molecule, removing each atom from the
    // board as it's discovered.
    while (board->movements.cursor < board->movements.length) {
        struct movement m = board->movements.elements[board->movements.cursor];
        for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
            if (!(m.atom & (BOND_LOW_BITS << bond_direction)))
                continue;
            struct vector p = m.position;
            struct vector d = v_offset_for_direction(bond_direction);
            p.u += d.u;
            p.v += d.v;
            atom *b = lookup_atom(board, p);
            if (!(*b & VALID) || (*b & REMOVED) || (*b & BEING_PRODUCED))
                continue;
            enqueue_movement(board, (struct movement){
                .position = p,
                .atom = *b,
                .base = base,
                .rotation = rotation,
                .translation = translation,
            });
            remove_atom(board, b);
        }
        board->movements.cursor++;
    }
}

static void perform_arm_instructions(struct solution *solution, struct board *board)
{
    board->movements.length = 0;
    board->movements.cursor = 0;
    uint32_t n = solution->number_of_arms;
    // this `ii` thing is a kind of cheesy way to do two passes; one pass for
    // drops and one for grabs (so handoffs work).
    for (uint32_t ii = 0; ii < 2*n; ++ii) {
        uint32_t i = ii;
        if (i >= n)
            i -= n;
        if (board->cycle < solution->arm_tape_start_cycle[i])
            continue;
        size_t index = board->cycle - solution->arm_tape_start_cycle[i];
        index %= solution->tape_period;
        if (index >= solution->arm_tape_length[i])
            continue;
        char inst = solution->arm_tape[i][index];
        if (inst == ' ' || inst == '\0')
            continue;
        // on half-cycle 1, grab and drop.
        // on half-cycle 2, perform the rest of the instructions.
        if ((board->half_cycle == 1) != (inst == 'r' || inst == 'f'))
            continue;
        // perform grabs after drops so handoffs work.
        if ((ii != i) != (inst == 'r'))
            continue;
        struct mechanism *m = &solution->arms[i];
        struct vector track_motion = {0, 0};
        // first, validate the instruction.
        if (inst == 'r') {
            // if the arm is already grabbing, grabbing again does nothing.
            // additionally, van berlo's wheel doesn't grab or release.
            if ((m->type & GRABBING) || (m->type & VAN_BERLO))
                continue;
        } else if (inst == 'f') {
            // if the arm isn't grabbing, then releasing does nothing.
            // additionally, van berlo's wheel doesn't grab or release.
            if (!(m->type & GRABBING) || (m->type & VAN_BERLO))
                continue;
        } else if ((inst == 'w' || inst == 's') && !(m->type & PISTON)) {
            // fprintf(stderr, "trying to extend/retract a non-piston arm");
            continue;
        } else if (inst == 'w') {
            // don't extend pistons past 3 hexes of length.
            if (m->direction_v.u == 3 || m->direction_v.u == -3 ||
             m->direction_v.v == 3 || m->direction_v.v == -3)
                continue;
        } else if (inst == 's') {
            // don't retract pistons below 1 hex of length.
            if (m->direction_v.u == 1 || m->direction_v.u == -1 ||
             m->direction_v.v == 1 || m->direction_v.v == -1)
                continue;
        } else if (inst == 't' || inst == 'g') {
            uint32_t index;
            if (!lookup_track(solution, m->position, &index)) {
                // fprintf(stderr, "trying to move an arm along a track that isn't on a track at %" PRId32 " %" PRId32 "\n", m->position.u, m->position.v);
                continue;
            }
            if (inst == 't')
                track_motion = solution->track_minus_motions[index];
            else if (inst == 'g')
                track_motion = solution->track_plus_motions[index];
            // if the motion amount is zero, this is the end of the track.
            if (track_motion.u == 0 && track_motion.v == 0)
                continue;
        }
        // next, apply the instruction to any grabbed atoms. atom movements are
        // added to a list (board->movements) and deferred until later.
        int step;
        switch (m->type & ANY_ARM) {
        case ARM:
        case PISTON:
            step = 6;
            break;
        case TWO_ARM:
            step = 3;
            break;
        case THREE_ARM:
            step = 2;
            break;
        case SIX_ARM:
        case VAN_BERLO:
        default:
            step = 1;
            break;
        }
        for (int direction = 0; direction < 6; direction += step) {
            struct vector offset = v_offset_for_direction(direction);
            atom *a = get_atom(board, *m, offset.u, offset.v);
            if (!*a)
                continue;
            if (inst == 'r') {
                atom grabs = (*a & GRABBED) / GRABBED_ONCE;
                *a &= ~GRABBED;
                *a |= (grabs + 1) * GRABBED_ONCE;
                m->type |= GRABBING_LOW_BIT << direction;
                continue;
            }
            if (!(m->type & (GRABBING_LOW_BIT << direction)) || !(*a & GRABBED))
                continue;
            if (inst == 'f') {
                atom grabs = (*a & GRABBED) / GRABBED_ONCE;
                *a &= ~GRABBED;
                *a |= (grabs - 1) * GRABBED_ONCE;
                continue;
            } 
            struct vector atom_pos = mechanism_relative_position(*m, offset.u, offset.v, 1);
            switch (inst) {
            case 'q': // pivot ccw
                move_atoms(board, a, atom_pos, atom_pos, 1, 0);
                break;
            case 'e': // pivot cw
                move_atoms(board, a, atom_pos, atom_pos, -1, 0);
                break;
            case 'a': // rotate ccw
                move_atoms(board, a, atom_pos, m->position, 1, 0);
                break;
            case 'd': // rotate cw
                move_atoms(board, a, atom_pos, m->position, -1, 0);
                break;
            case 'w': // extend piston
                move_atoms(board, a, atom_pos, normalize_axis(m->direction_v), 0, 1);
                break;
            case 's': // retract piston
                move_atoms(board, a, atom_pos, normalize_axis(m->direction_v), 0, -1);
                break;
            case 't': // move along track, - direction
            case 'g': // move along track, + direction
                move_atoms(board, a, atom_pos, track_motion, 0, 1);
                break;
            default:
                break;
            }
        }
        // finally, transform the arm itself.
        switch (inst) {
        case 'a': { // rotate ccw
            struct vector u = u_offset_for_direction(1);
            struct vector v = v_offset_for_direction(1);
            struct vector nu = {
                m->direction_u.u * u.u + m->direction_v.u * u.v,
                m->direction_u.v * u.u + m->direction_v.v * u.v,
            };
            struct vector nv = {
                m->direction_u.u * v.u + m->direction_v.u * v.v,
                m->direction_u.v * v.u + m->direction_v.v * v.v,
            };
            m->direction_u = nu;
            m->direction_v = nv;
            break;
        }
        case 'd': { // rotate cw
            struct vector u = u_offset_for_direction(-1);
            struct vector v = v_offset_for_direction(-1);
            struct vector nu = {
                m->direction_u.u * u.u + m->direction_v.u * u.v,
                m->direction_u.v * u.u + m->direction_v.v * u.v,
            };
            struct vector nv = {
                m->direction_u.u * v.u + m->direction_v.u * v.v,
                m->direction_u.v * v.u + m->direction_v.v * v.v,
            };
            m->direction_u = nu;
            m->direction_v = nv;
            break;
        }
        case 'w': // extend piston
            adjust_axis_magnitude(&m->direction_u, 1);
            adjust_axis_magnitude(&m->direction_v, 1);
            break;
        case 's': // retract piston
            adjust_axis_magnitude(&m->direction_u, -1);
            adjust_axis_magnitude(&m->direction_v, -1);
            break;
        case 't': // move along track, - direction
        case 'g': // move along track, + direction
            m->position.u += track_motion.u;
            m->position.v += track_motion.v;
            break;
        case 'r': // grab
            m->type |= GRABBING;
            break;
        case 'f': // drop
            // clear all the mechanism grabbing bits at once.
            m->type &= ~GRABBING_EVERYTHING;
            m->type &= ~GRABBING;
            break;
        default:
            break;
        }
    }
    // carry out deferred atom movements.
    ensure_capacity(board, board->movements.length);
    for (size_t i = 0; i < board->movements.length; ++i) {
        struct movement m = board->movements.elements[i];
        // printf("movement: %d %d by %d / %d around %d %d\n", m.position.u, m.position.v, m.rotation, m.translation, m.base.u, m.base.v);
        struct vector delta = m.position;
        delta.u -= m.base.u;
        delta.v -= m.base.v;
        struct vector u = u_offset_for_direction(m.rotation);
        struct vector v = v_offset_for_direction(m.rotation);
        struct vector to = {
            u.u * delta.u + v.u * delta.v + m.base.u + m.translation * m.base.u,
            u.v * delta.u + v.v * delta.v + m.base.v + m.translation * m.base.v,
        };
        // rotate the atom's bonds by rotating the bits which represent those
        // bonds.
        atom bonds = m.atom & ALL_BONDS;
        // first, do a standard bit shift by the rotation amount.
        m.rotation = normalize_direction(m.rotation);
        bonds <<= m.rotation;
        // take the overflow bits that were shifted off the end...
        atom mask = 0x3FULL >> (6 - m.rotation);
        atom overflow = bonds & ((mask * BOND_LOW_BITS) << 6);
        bonds &= ~overflow;
        // ...and shift them back around to the other side.
        bonds |= overflow >> 6;
        m.atom &= ~ALL_BONDS;
        m.atom |= bonds;
        atom *a = insert_atom(board, to, "atom moved on top of another atom");
        *a = VALID | m.atom;

        // xx collision detection / error handling
    }
    if (board->half_cycle == 2) {
        for (uint32_t i = 0; i < n; ++i) {
            struct mechanism m = solution->arms[i];
            atom a = *lookup_atom(board, m.position);
            if ((a & VALID) && !(a & REMOVED))
                report_collision(board, m.position, "arm base collided with atom");
        }
    }
}

static void flag_blocked_inputs(struct solution *solution, struct board *board)
{
    uint32_t m = solution->number_of_input_molecules;
    for (uint32_t i = 0; i < m; ++i)
        solution->input_molecule_is_blocked[i] = false;
    for (uint32_t i = 0, molecule = 0; i < m; ++i) {
        uint32_t n = solution->input_molecule_lengths[i];
        for (uint32_t j = 0; j < n; ++j) {
            struct vector p = solution->input_atom_positions[molecule + j];
            atom a = *lookup_atom(board, p);
            if ((a & VALID) && !(a & REMOVED) && !(a & BEING_PRODUCED)) {
                solution->input_molecule_is_blocked[i] = true;
                break;
            }
        }
        molecule += n;
    }
}

static void spawn_inputs(struct solution *solution, struct board *board)
{
    uint32_t m = solution->number_of_input_molecules;
    for (uint32_t i = 0, molecule = 0; i < m; ++i) {
        uint32_t n = solution->input_molecule_lengths[i];
        if (!solution->input_molecule_is_blocked[i]) {
            ensure_capacity(board, n);
            for (uint32_t j = 0; j < n; ++j) {
                atom input = solution->input_atoms[molecule + j];
                struct vector p = solution->input_atom_positions[molecule + j];
                atom *a = insert_atom(board, p, "overlapped inputs");
                *a = VALID | input;
            }
        }
        molecule += n;
    }
}

static int output_value;
static bool wrong_output;

static void consume_outputs(struct solution *solution, struct board *board)
{
    uint32_t m = solution->number_of_output_molecules;
    for (uint32_t i = 0, molecule = 0; i < m; ++i) {
        uint32_t n = solution->output_molecule_lengths[i];
        // first, check the entire output to see if it matches.
        bool match = true;
        wrong_output = true;
        output_value = 0;
        for (uint32_t j = 0; j < n; ++j) {
            struct vector p = solution->output_atom_positions[molecule + j];
            atom output = solution->output_atoms[molecule + j];
            atom *a = lookup_atom(board, p);
            if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED)) {
                match = false;
                wrong_output = false;
                break;
            }
            if ((*a & ALL_BONDS) != (output & ALL_BONDS) || (*a & GRABBED)) {
                match = false;
                wrong_output = false;
                break;
            }
            int k = j;
            if (j == 4)
                k = 5;
            else if (j == 5)
                k = 4;
            output_value <<= 1;
            output_value |= !(*a & SALT);
            if ((*a & (ANY_ATOM | ALL_BONDS)) != output) {
                // this could be a "wrong output".
                match = false;
            }
        }
        // if the output is a match, remove it and increment the output counter.
        if (match) {
            wrong_output = false;
            for (uint32_t j = 0; j < n; ++j) {
                struct vector p = solution->output_atom_positions[molecule + j];
                atom *a = lookup_atom(board, p);
                remove_atom(board, a);
            }
            solution->output_molecule_number_of_outputs[i]++;
        } else if (wrong_output) {
            board->collision = true;
            board->collision_location = solution->output_atom_positions[molecule];
            board->collision_reason = "output didn't match";
            break;
        }
        molecule += n;
    }
}

static void flag_unbonded_atoms(struct board *board)
{
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom *a = &board->atoms_at_positions[i].atom;
        if (!(*a & VALID) || (*a & REMOVED))
            continue;
        if (!(*a & ALL_BONDS))
            *a |= UNBONDED;
        else
            *a &= ~UNBONDED;
    }
}

static bool check_completion(struct solution *solution)
{
    uint64_t min = UINT64_MAX;
    uint32_t m = solution->number_of_output_molecules;
    for (uint32_t i = 0; i < m; ++i) {
        uint64_t count = solution->output_molecule_number_of_outputs[i];
        if (count < min)
            min = count;
    }
    return min >= solution->target_number_of_outputs;
}

static int _last_cycle = 0;
static int _last_capacity = 0;
static int _last_query_u = 0;
static int _last_query_v = 0;
static int _last_used = 0;
static int _last_removed = 0;
int last_cycle(void)
{
    return _last_cycle;
}
int last_capacity(void)
{
    return _last_capacity;
}
int last_query_u(void)
{
    return _last_query_u;
}
int last_query_v(void)
{
    return _last_query_v;
}
int last_used(void)
{
    return _last_used;
}
int last_removed(void)
{
    return _last_removed;
}

static void cycle(struct solution *solution, struct board *board)
{
    for (board->half_cycle = 1; board->half_cycle <= 2; board->half_cycle++) {
        perform_arm_instructions(solution, board);
        flag_blocked_inputs(solution, board);
        spawn_inputs(solution, board);
        flag_unbonded_atoms(board);
        apply_glyphs(solution, board);
        consume_outputs(solution, board);
    }
    board->complete = check_completion(solution);
    board->cycle++;
    _last_cycle = board->cycle;
}

static void create_van_berlo_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv, atom element)
{
    ensure_capacity(board, 1);
    struct vector p = mechanism_relative_position(m, du, dv, 1);
    atom *a = insert_atom(board, p, "van berlo overlap");
    *a = VALID | GRABBED_ONCE | VAN_BERLO | element;
}

#if 0
static void print_board(struct board *board)
{
    for (uint32_t i = 0; i < board->capacity; ++i) {
        atom a = board->atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector position = board->atoms_at_positions[i].position;
        printf("%" PRId32 " %" PRId32 " %" PRIx64 "\n", position.u, position.v, a);
    }
}
#endif

// see appendix for solution decoding.
static bool decode_solution(struct solution *solution, struct puzzle_file *pf,
 struct solution_file *sf);

static enum mechanism_type decode_mechanism_type(struct byte_string part_name);

// wk9 specific
static struct puzzle_file *pf;
static struct solution_file *sf;
static int test_one(int a, int b, bool fun)
{
    struct solution solution = { 0 };
    struct board board = { 0 };

    if (!decode_solution(&solution, pf, sf))
        return -999;

    if (solution.number_of_output_molecules != 1) {
        fprintf(stderr, "output must be placed\n");
        return -999;
    }

    uint32_t atom_index = 0;
    for (uint32_t i = 0; i < solution.number_of_input_molecules; ++i) {
        if (solution.input_molecule_puzzle_index[i] == 0) {
            for (int j = 0; j < 8; ++j) {
                solution.input_atoms[atom_index + 7 - j] &= ~ANY_ATOM;
                solution.input_atoms[atom_index + 7 - j] |= ((a >> j) & 1) ? FIRE : SALT;
            }
        } else if (solution.input_molecule_puzzle_index[i] == 1) {
            for (int j = 0; j < 8; ++j) {
                solution.input_atoms[atom_index + 7 - j] &= ~ANY_ATOM;
                solution.input_atoms[atom_index + 7 - j] |= ((b >> j) & 1) ? FIRE : SALT;
            }
        }
        atom_index += solution.input_molecule_lengths[i];
    }
    for (int i = 0; i < 8; ++i) {
        int j = i;
        if (i == 4)
            j = 5;
        else if (i == 5)
            j = 4;
        solution.output_atoms[7 - j] &= ~ANY_ATOM;
        solution.output_atoms[7 - j] |= (((a - b) >> i) & 1) ? FIRE : SALT;
    }

    ensure_capacity(&board, 1);
    spawn_inputs(&solution, &board);
    for (size_t i = 0; i < solution.number_of_arms; ++i) {
        if (!(solution.arms[i].type & VAN_BERLO))
            continue;
        solution.arms[i].type |= GRABBING;
        solution.arms[i].type |= GRABBING_EVERYTHING;
        create_van_berlo_atom(&board, solution.arms[i], 0, 1, SALT);
        create_van_berlo_atom(&board, solution.arms[i], 1, 0, WATER);
        create_van_berlo_atom(&board, solution.arms[i], 1, -1, AIR);
        create_van_berlo_atom(&board, solution.arms[i], 0, -1, SALT);
        create_van_berlo_atom(&board, solution.arms[i], -1, 0, FIRE);
        create_van_berlo_atom(&board, solution.arms[i], -1, 1, EARTH);
    }

    while (board.cycle < 10000 * (1 + solution.output_molecule_number_of_outputs[0]) && !board.complete) {
        cycle(&solution, &board);
        if (board.collision) {
            if (!fun) {
                fprintf(stderr, "collision on cycle %llu at %" PRId32 ", %" PRId32 ": %s\n", board.cycle,
                 board.collision_location.u, board.collision_location.v,
                 board.collision_reason);
            }
            break;
        }
    }
    fflush(stdout);
    fflush(stderr);

    free(solution.glyphs);
    free(solution.arms);
    for (uint32_t i = 0; i < solution.number_of_arms; ++i)
        free(solution.arm_tape[i]);
    free(solution.arm_tape);
    free(solution.arm_tape_length);
    free(solution.arm_tape_start_cycle);
    free(solution.track_positions);
    free(solution.track_plus_motions);
    free(solution.track_minus_motions);
    free(solution.input_atoms);
    free(solution.input_atom_positions);
    free(solution.input_molecule_lengths);
    free(solution.input_molecule_is_blocked);
    free(solution.input_molecule_puzzle_index);
    free(solution.output_atoms);
    free(solution.output_atom_positions);
    free(solution.output_molecule_lengths);
    free(solution.output_molecule_number_of_outputs);

    free(board.atoms_at_positions);
    free(board.movements.elements);

    if (wrong_output) {
        // fix output value.
        int bit4 = (output_value >> 4) & 1;
        int bit5 = (output_value >> 5) & 1;
        int fixed = output_value;
        fixed &= ~(0x3 << 4);
        fixed |= bit4 << 5;
        fixed |= bit5 << 4;
        return -fixed;
    }
    if (!board.complete)
        return -998;
    return board.cycle;
    // printf("solution file says cycle count is: %" PRIu32 "\n", sf->cycles);
    // printf("simulation says cycle count is: %" PRIu64 "\n", board.cycle);
}

int test(int a, int b, int fun)
{
    if (fun) {
        double x = 0;
        double y = 0;
        uint32_t n = 0;
        for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
            struct byte_string saved = sf->parts[i].name;
            enum mechanism_type type = decode_mechanism_type(saved);
            if (type == 0 || type == EQUILIBRIUM)
                continue;
            n++;
        }
        uint32_t j = 0;
        for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
            struct byte_string saved = sf->parts[i].name;
            enum mechanism_type type = decode_mechanism_type(saved);
            if (type == 0 || type == EQUILIBRIUM)
                continue;
            sf->parts[i].name = (struct byte_string){ (unsigned char *)"disabled", strlen("disabled") };
            int result = test_one(a, b, true);
            if (result > 0) {
                x += cos(M_PI * 2 * (double)j / (double)n);
                y += sin(M_PI * 2 * (double)j / (double)n);
            }
            j++;
            sf->parts[i].name = saved;
        }
        if (x == 0 && y == 0)
            return -998;
        return -1000 - (int)(65536 * (atan2(y, x) / M_PI + 1) / 2);
    } else
        return test_one(a, b, false);
}

static unsigned char OM2021_W9_00_puzzle_bytes[] = {
  0x03, 0x00, 0x00, 0x00, 0x14, 0x45, 0x58, 0x50, 0x4c, 0x4f, 0x53, 0x49,
  0x56, 0x45, 0x20, 0x4c, 0x4f, 0x47, 0x49, 0x43, 0x20, 0x55, 0x4e, 0x49,
  0x54, 0x62, 0x7f, 0x8b, 0x02, 0x01, 0x00, 0x10, 0x01, 0x0f, 0x3f, 0xc7,
  0x17, 0x00, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00,
  0x00, 0x01, 0xfd, 0x00, 0x01, 0xfe, 0x00, 0x01, 0xff, 0x00, 0x04, 0x00,
  0x00, 0x01, 0x01, 0x00, 0x04, 0x02, 0x00, 0x04, 0x03, 0x00, 0x01, 0x04,
  0x00, 0x07, 0x00, 0x00, 0x00, 0x01, 0xfd, 0x00, 0xfe, 0x00, 0x01, 0xfe,
  0x00, 0xff, 0x00, 0x01, 0xff, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x01,
  0x00, 0x01, 0x01, 0x00, 0x02, 0x00, 0x01, 0x02, 0x00, 0x03, 0x00, 0x01,
  0x03, 0x00, 0x04, 0x00, 0x08, 0x00, 0x00, 0x00, 0x01, 0xfd, 0x00, 0x01,
  0xfe, 0x00, 0x01, 0xff, 0x00, 0x04, 0x00, 0x00, 0x01, 0x01, 0x00, 0x01,
  0x02, 0x00, 0x04, 0x03, 0x00, 0x04, 0x04, 0x00, 0x07, 0x00, 0x00, 0x00,
  0x01, 0xfd, 0x00, 0xfe, 0x00, 0x01, 0xfe, 0x00, 0xff, 0x00, 0x01, 0xff,
  0x00, 0x00, 0x00, 0x01, 0x02, 0x00, 0x03, 0x00, 0x01, 0x03, 0x00, 0x04,
  0x00, 0x01, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x02,
  0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x01, 0x01, 0x00, 0x01, 0x00, 0x00,
  0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x02, 0x00, 0x00, 0x00, 0x04, 0x00,
  0x00, 0x01, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x01,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, 0x01, 0xfd, 0x00,
  0x01, 0xfe, 0x00, 0x01, 0x00, 0x00, 0x01, 0xff, 0x00, 0x01, 0x01, 0x00,
  0x01, 0x02, 0x00, 0x04, 0x03, 0x00, 0x04, 0x04, 0x00, 0x07, 0x00, 0x00,
  0x00, 0x01, 0x03, 0x00, 0x04, 0x00, 0x01, 0x02, 0x00, 0x03, 0x00, 0x01,
  0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x01, 0xff, 0x00,
  0x00, 0x00, 0x01, 0xfe, 0x00, 0xff, 0x00, 0x01, 0xfd, 0x00, 0xfe, 0x00,
  0x01, 0x00, 0x00, 0x00, 0x00
};
static struct byte_string OM2021_W9_00_puzzle =
    { OM2021_W9_00_puzzle_bytes, sizeof(OM2021_W9_00_puzzle_bytes) };

void init(void)
{
    pf = parse_puzzle_byte_string(OM2021_W9_00_puzzle);
}

void load_solution(unsigned char *bytes, size_t len)
{
    if (sf)
        free_solution_file(sf);
    struct byte_string b = { malloc(len), len };
    memcpy(b.bytes, bytes, len);
    sf = parse_solution_byte_string(b);
    if (!sf) {
        fprintf(stderr, "not a valid solution file\n");
        free(b.bytes);
    }
    sf->owns_bytes = true;
}

int main(int argc, char *argv[])
{
    // read input files.
    if (argc < 2) {
        fprintf(stderr, "usage: wk9 [solution]\n");
        return -1;
    }
    init();
    sf = parse_solution_file(argv[1]);
    if (!sf) {
        fprintf(stderr, "%s is not a valid solution file\n", argv[1]);
        return -1;
    }
    for (int a = -128; a < 127; ++a) {
        for (int b = -128; b < 127; ++b) {
            if (a - b < -128 || a - b > 127)
                continue;
            printf("%d %d %d\n", a, b, test(a, b, 1));
        }
    }
}

// appendix A -- hash table functions.

// This is the 32-bit FNV-1a hash diffusion algorithm.
// http://www.isthe.com/chongo/tech/comp/fnv/index.html

static uint32_t fnv(const void *dataPointer, size_t length)
{
    uint32_t hash = 0x811c9dc5;
    const unsigned char *data = dataPointer;
    for (size_t i = 0; i < length; ++i) {
        hash ^= data[i];
        hash *= 0x01000193;
    }
    return hash;
}

static bool lookup_track(struct solution *solution, struct vector query, uint32_t *index)
{
    if (solution->track_table_size == 0)
        return false;
    uint32_t hash = fnv(&query, sizeof(query));
    uint32_t mask = solution->track_table_size - 1;
    *index = hash & mask;
    while (true) {
        struct vector p = solution->track_positions[*index];
        if (p.u == INT32_MIN && p.v == INT32_MIN)
            return false;
        if (p.u == query.u && p.v == query.v)
            return true;
        *index = (*index + 1) & mask;
        if (*index == (hash & mask))
            abort();
    }
}

static atom *lookup_atom(struct board *board, struct vector query)
{
    _last_used = board->used;
    _last_removed = board->removed;
    _last_capacity = board->capacity;
    _last_query_u = query.u;
    _last_query_v = query.v;
    uint32_t hash = fnv(&query, sizeof(query));
    uint32_t mask = board->capacity - 1;
    uint32_t index = hash & mask;
    uint32_t steps = 0;
    while (true) {
        atom *a = &board->atoms_at_positions[index].atom;
        if (!(*a & VALID))
            return a;
        struct vector position = board->atoms_at_positions[index].position;
        if (!memcmp(&position, &query, sizeof(position)))
            return a;
        index = (index + 1) & mask;
        steps++;
        if (index == (hash & mask))
            abort();
    }
}

static atom *insert_atom(struct board *board, struct vector query, const char *collision_reason)
{
    _last_used = board->used;
    _last_removed = board->removed;
    _last_capacity = board->capacity;
    _last_query_u = query.u;
    _last_query_v = query.v;
    uint32_t hash = fnv(&query, sizeof(query));
    uint32_t mask = board->capacity - 1;
    uint32_t index = hash & mask;
    uint32_t removed = UINT32_MAX;
    while (true) {
        atom *a = &board->atoms_at_positions[index].atom;
        // in order to detect overlapping atoms, continue searching even if
        // there's a removed atom in the table.
        if (removed == UINT32_MAX && (*a & REMOVED))
            removed = index;
        if (!(*a & VALID)) {
            if (removed != UINT32_MAX) {
                board->removed--;
                index = removed;
            }
            board->atoms_at_positions[index].position = query;
            break;
        }
        struct vector position = board->atoms_at_positions[index].position;
        if (position.u == query.u && position.v == query.v && !(*a & REMOVED))
            break;
        index = (index + 1) & mask;
        if (index == (hash & mask)) {
            if (removed != UINT32_MAX) {
                index = removed;
                board->removed--;
                board->atoms_at_positions[index].position = query;
                break;
            }
            abort();
        }
    }
    atom *a = &board->atoms_at_positions[index].atom;
    if ((*a & VALID) && !(*a & REMOVED))
        report_collision(board, query, collision_reason);
    board->used++;
    return a;
}

static void ensure_capacity(struct board *board, uint32_t potential_insertions)
{
    uint32_t n = board->capacity;
    if (n == 0)
        n = 16;
    while (2 * n <= 7 * (board->used + potential_insertions))
        n *= 2;
    if (n == board->capacity && 3 * (board->removed + board->used) < 2 * board->capacity)
        return;
    struct board old = *board;
    board->atoms_at_positions = calloc(n, sizeof(*board->atoms_at_positions));
    board->capacity = n;
    board->used = 0;
    board->removed = 0;
    for (uint32_t i = 0; i < old.capacity; ++i) {
        atom a = old.atoms_at_positions[i].atom;
        if (!(a & VALID) || (a & REMOVED))
            continue;
        struct vector position = old.atoms_at_positions[i].position;
        *insert_atom(board, position, "reinserting colliding atoms") = a;
    }
    free(old.atoms_at_positions);
}

// appendix B -- puzzle / solution file decoding.

static char decode_instruction(char inst)
{
    switch (inst) {
    case 'R': return 'd';
    case 'r': return 'a';
    case 'E': return 'w';
    case 'e': return 's';
    case 'G': return 'r';
    case 'g': return 'f';
    case 'P': return 'e';
    case 'p': return 'q';
    case 'A': return 'g';
    case 'a': return 't';
    case 'C': fprintf(stderr, "repeat instruction unsupported\n"); return ' ';
    case 'X': fprintf(stderr, "reset instruction unsupported\n"); return ' ';
    case 'O': return ' ';
    case ' ': return ' ';
    default: return ' ';
    }
}

static atom decode_atom(uint32_t atom)
{
    return 1ULL << atom;
}

static atom decode_bond_type(unsigned char bond_type)
{
    atom bonds = 0;
    if (bond_type & 1)
        bonds |= NORMAL_BONDS;
    if (bond_type & 2)
        bonds |= TRIPLEX_R_BONDS;
    if (bond_type & 4)
        bonds |= TRIPLEX_K_BONDS;
    if (bond_type & 8)
        bonds |= TRIPLEX_Y_BONDS;
    return bonds;
}

static enum mechanism_type decode_mechanism_type(struct byte_string part_name)
{
    if (byte_string_is(part_name, "glyph-calcification"))
        return CALCIFICATION;
    else if (byte_string_is(part_name, "glyph-life-and-death"))
        return ANIMISMUS;
    else if (byte_string_is(part_name, "glyph-projection"))
        return PROJECTION;
    else if (byte_string_is(part_name, "glyph-dispersion"))
        return DISPERSION;
    else if (byte_string_is(part_name, "glyph-purification"))
        return PURIFICATION;
    else if (byte_string_is(part_name, "glyph-duplication"))
        return DUPLICATION;
    else if (byte_string_is(part_name, "glyph-unification"))
        return UNIFICATION;
    else if (byte_string_is(part_name, "bonder"))
        return BONDING;
    else if (byte_string_is(part_name, "unbonder"))
        return UNBONDING;
    else if (byte_string_is(part_name, "bonder-prisma"))
        return TRIPLEX_BONDING;
    else if (byte_string_is(part_name, "bonder-speed"))
        return MULTI_BONDING;
    else if (byte_string_is(part_name, "glyph-disposal"))
        return DISPOSAL;
    else if (byte_string_is(part_name, "glyph-marker"))
        return EQUILIBRIUM;
    else if (byte_string_is(part_name, "arm1"))
        return ARM;
    else if (byte_string_is(part_name, "arm2"))
        return TWO_ARM;
    else if (byte_string_is(part_name, "arm3"))
        return THREE_ARM;
    else if (byte_string_is(part_name, "arm6"))
        return SIX_ARM;
    else if (byte_string_is(part_name, "piston"))
        return PISTON;
    else if (byte_string_is(part_name, "baron"))
        return VAN_BERLO;
    else
        return 0;
}

static void decode_molecule(struct puzzle_molecule c, struct mechanism m, atom *atoms, struct vector *positions)
{
    for (uint32_t i = 0; i < c.number_of_atoms; ++i) {
        atoms[i] = decode_atom(c.atoms[i].type);
        positions[i] = mechanism_relative_position(m, c.atoms[i].offset[1], c.atoms[i].offset[0], 1);
    }
    for (uint32_t i = 0; i < c.number_of_bonds; ++i) {
        struct puzzle_bond b = c.bonds[i];
        atom bonds = decode_bond_type(c.bonds[i].type);
        struct vector p1 = mechanism_relative_position(m, b.from[1], b.from[0], 1);
        struct vector p2 = mechanism_relative_position(m, b.to[1], b.to[0], 1);
        atom b1 = bonds & bond_direction(m, b.to[1] - b.from[1], b.to[0] - b.from[0]);
        atom b2 = bonds & bond_direction(m, b.from[1] - b.to[1], b.from[0] - b.to[0]);
        // yes, this is O(n^2).
        for (uint32_t j = 0; j < c.number_of_atoms; ++j) {
            struct vector p = positions[j];
            if (p.u == p1.u && p.v == p1.v)
                atoms[j] |= b1;
            if (p.u == p2.u && p.v == p2.v)
                atoms[j] |= b2;
        }
    }
}

static bool decode_solution(struct solution *solution, struct puzzle_file *pf, struct solution_file *sf)
{
    size_t number_of_track_hexes = 0;
    uint32_t number_of_input_atoms = 0;
    uint32_t number_of_output_atoms = 0;
    // first pass through the solution file: count how many things of each type
    // there are.  these counts are used to allocate arrays of the correct size.
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        enum mechanism_type type = decode_mechanism_type(sf->parts[i].name);
        if (type & ANY_ARM)
            solution->number_of_arms++;
        else if (type & ANY_GLYPH)
            solution->number_of_glyphs++;
        else if (byte_string_is(sf->parts[i].name, "track"))
            number_of_track_hexes += sf->parts[i].number_of_track_hexes;
        else if (byte_string_is(sf->parts[i].name, "input")) {
            if (sf->parts[i].which_input_or_output >= pf->number_of_inputs) {
                fprintf(stderr, "solution file error: input out of range\n");
                return false;
            }
            solution->number_of_input_molecules++;
            number_of_input_atoms += pf->inputs[sf->parts[i].which_input_or_output].number_of_atoms;
        } else if (byte_string_is(sf->parts[i].name, "out-std")) {
            if (sf->parts[i].which_input_or_output >= pf->number_of_outputs) {
                fprintf(stderr, "solution file error: output out of range\n");
                return false;
            }
            solution->number_of_output_molecules++;
            number_of_output_atoms += pf->outputs[sf->parts[i].which_input_or_output].number_of_atoms;
        } else if (byte_string_is(sf->parts[i].name, "out-rep")) {
            // todo
            fprintf(stderr, "solution file error: repeating outputs are currently unsupported\n");
            return false;
        } else if (byte_string_is(sf->parts[i].name, "pipe")) {
            // todo
            fprintf(stderr, "solution file error: conduits are currently unsupported\n");
            return false;
        }
    }
    // now that we know how many elements each array should have, allocate them
    // all here.
    solution->glyphs = calloc(solution->number_of_glyphs, sizeof(struct mechanism));

    solution->arms = calloc(solution->number_of_arms, sizeof(struct mechanism));
    solution->arm_tape = calloc(solution->number_of_arms, sizeof(char *));
    solution->arm_tape_length = calloc(solution->number_of_arms, sizeof(size_t));
    solution->arm_tape_start_cycle = calloc(solution->number_of_arms, sizeof(uint64_t));

    solution->input_atoms = calloc(number_of_input_atoms, sizeof(atom));
    solution->input_atom_positions = calloc(number_of_input_atoms, sizeof(struct vector));
    solution->input_molecule_lengths = calloc(solution->number_of_input_molecules, sizeof(uint32_t));
    solution->input_molecule_is_blocked = calloc(solution->number_of_input_molecules, sizeof(bool));
    solution->input_molecule_puzzle_index = calloc(solution->number_of_input_molecules, sizeof(uint32_t));

    solution->output_atoms = calloc(number_of_output_atoms, sizeof(atom));
    solution->output_atom_positions = calloc(number_of_output_atoms, sizeof(struct vector));
    solution->output_molecule_lengths = calloc(solution->number_of_output_molecules, sizeof(uint32_t));
    solution->output_molecule_number_of_outputs = calloc(solution->number_of_output_molecules, sizeof(uint64_t));

    solution->track_table_size = 1;
    while (2 * solution->track_table_size < 3 * number_of_track_hexes)
        solution->track_table_size *= 2;
    solution->track_positions = calloc(solution->track_table_size, sizeof(solution->track_positions[0]));
    for (int i = 0; i < solution->track_table_size; ++i)
        solution->track_positions[i] = (struct vector){ INT32_MIN, INT32_MIN };
    solution->track_plus_motions = calloc(solution->track_table_size, sizeof(solution->track_plus_motions[0]));
    solution->track_minus_motions = calloc(solution->track_table_size, sizeof(solution->track_minus_motions[0]));

    // second pass: fill in the arrays with the data from the file.
    uint32_t arm_index = 0;
    uint32_t glyph_index = 0;
    size_t input_molecule_index = 0;
    size_t input_atom_index = 0;
    size_t output_molecule_index = 0;
    size_t output_atom_index = 0;
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        struct solution_part part = sf->parts[i];
        struct mechanism m = {
            .type = decode_mechanism_type(part.name),
            .position = { part.position[1], part.position[0] },
            .direction_u = u_offset_for_direction(part.rotation),
            .direction_v = v_offset_for_direction(part.rotation),
        };
        m.direction_u.u *= part.size;
        m.direction_u.v *= part.size;
        m.direction_v.u *= part.size;
        m.direction_v.v *= part.size;
        if (m.type & ANY_ARM)
            solution->arms[arm_index++] = m;
        else if (m.type & ANY_GLYPH)
            solution->glyphs[glyph_index++] = m;
        else if (byte_string_is(part.name, "track")) {
            struct vector last_position = m.position;
            for (uint32_t j = 0; j < part.number_of_track_hexes + 1; ++j) {
                struct solution_hex_offset hex;
                if (j < part.number_of_track_hexes)
                    hex = part.track_hexes[j];
                else {
                    // two-hex tracks can't become a loop.
                    if (part.number_of_track_hexes <= 2)
                        break;
                    hex = part.track_hexes[0];
                    int32_t du = hex.offset[1] - part.track_hexes[j - 1].offset[1];
                    int32_t dv = hex.offset[0] - part.track_hexes[j - 1].offset[0];
                    // if the offset between the two hexes isn't a cardinal
                    // direction, then the ends are too far away for the track
                    // to become a loop.
                    if (direction_for_offset((struct vector){ du, dv }) < 0)
                        break;
                }
                struct vector p = mechanism_relative_position(m, hex.offset[1], hex.offset[0], 1);
                uint32_t index;
                lookup_track(solution, p, &index);
                solution->track_positions[index] = p;
                if (j != 0) {
                    solution->track_minus_motions[index] = (struct vector){ last_position.u - p.u, last_position.v - p.v };
                    lookup_track(solution, last_position, &index);
                    solution->track_plus_motions[index] = (struct vector){ p.u - last_position.u, p.v - last_position.v };
                }
                last_position = p;
            }
        } else if (byte_string_is(part.name, "input")) {
            struct puzzle_molecule c = pf->inputs[part.which_input_or_output];
            solution->input_molecule_lengths[input_molecule_index] = c.number_of_atoms;
            decode_molecule(c, m, solution->input_atoms + input_atom_index,
             solution->input_atom_positions + input_atom_index);
            solution->input_molecule_puzzle_index[input_molecule_index] = part.which_input_or_output;
            input_atom_index += c.number_of_atoms;
            input_molecule_index++;
        } else if (byte_string_is(part.name, "out-std")) {
            struct puzzle_molecule c = pf->outputs[part.which_input_or_output];
            solution->output_molecule_lengths[output_molecule_index] = c.number_of_atoms;
            decode_molecule(c, m, solution->output_atoms + output_atom_index,
             solution->output_atom_positions + output_atom_index);
            output_atom_index += c.number_of_atoms;
            output_molecule_index++;
        }
    }
    // decode arm tapes in one final pass.  this has to be another pass because
    // reset instructions depend on where track has been placed.
    arm_index = 0;
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        struct solution_part part = sf->parts[i];
        if (!(decode_mechanism_type(part.name) & ANY_ARM))
            continue;
        uint32_t min_tape = UINT32_MAX;
        uint32_t max_tape = 0;
        for (uint32_t j = 0; j < part.number_of_instructions; ++j) {
            uint32_t index = part.instructions[j].index;
            if (index < min_tape)
                min_tape = index;
            if (index > max_tape)
                max_tape = index;
        }
        if (max_tape < min_tape) {
            arm_index++;
            continue;
        }
        uint32_t tape_length = max_tape - min_tape + 1;
        // multiply by two to leave room for reset and repeat instructions
        // (which cannot exceed the length of the earlier instructions even in
        // the worst case).
        solution->arm_tape[arm_index] = calloc(tape_length * 2, 1);
        solution->arm_tape_start_cycle[arm_index] = min_tape;
        uint32_t last_end = 0;
        uint32_t last_repeat = 0;
        uint32_t reset_from = 0;
        char *tape = solution->arm_tape[arm_index];
        for (uint32_t j = 0; j < part.number_of_instructions; ++j) {
            struct solution_instruction inst = part.instructions[j];
            uint32_t n = inst.index - min_tape;
            if (inst.instruction == 'C') { // repeat
                while (j < part.number_of_instructions && part.instructions[j].instruction == 'C') {
                    memcpy(tape + part.instructions[j].index - min_tape, tape + last_repeat, last_end - last_repeat);
                    uint32_t m = part.instructions[j].index - min_tape + last_end - last_repeat;
                    if (m > tape_length)
                        tape_length = m;
                    j++;
                }
                if (j < part.number_of_instructions) {
                    last_repeat = part.instructions[j].index - min_tape;
                    reset_from = last_repeat;
                }
                j--;
            } else if (inst.instruction == 'X') { // reset
                struct vector position = solution->arms[arm_index].position;
                struct vector track = position;
                int track_steps = 0;
                int rotation = 0;
                int piston = part.size;
                int grab = 0;
                for (uint32_t k = reset_from; k < n; ++k) {
                    struct vector motion = { 0, 0 };
                    if (tape[k] == 'a') {
                        rotation++;
                        if (rotation > 3)
                            rotation -= 6;
                    } else if (tape[k] == 'd') {
                        rotation--;
                        if (rotation < -3)
                            rotation += 6;
                    } else if (tape[k] == 'w') {
                        if (piston < 3)
                            piston++;
                    } else if (tape[k] == 's') {
                        if (piston > 1)
                            piston--;
                    } else if (tape[k] == 'g') {
                        uint32_t index;
                        if (lookup_track(solution, track, &index)) {
                            motion = solution->track_plus_motions[index];
                            if (motion.u != 0 || motion.v != 0)
                                track_steps++;
                        }
                    } else if (tape[k] == 't') {
                        uint32_t index;
                        if (lookup_track(solution, track, &index)) {
                            motion = solution->track_minus_motions[index];
                            if (motion.u != 0 || motion.v != 0)
                                track_steps--;
                        }
                    }
                    else if (tape[k] == 'r')
                        grab = 1;
                    else if (tape[k] == 'f')
                        grab = 0;
                    track.u += motion.u;
                    track.v += motion.v;
                    if (track.u == position.u && track.v == position.v)
                        track_steps = 0;
                }
                if (grab > 0)
                    tape[n++] = 'f';
                while (piston > part.size) {
                    tape[n++] = 's';
                    piston--;
                }
                while (piston < part.size) {
                    tape[n++] = 'w';
                    piston++;
                }
                while (rotation > 0) {
                    tape[n++] = 'd';
                    rotation--;
                }
                while (rotation < 0) {
                    tape[n++] = 'a';
                    rotation++;
                }
                if (track_steps != 0) {
                    // look for a path forward on the track that's shorter than
                    // the path backward.
                    int search_depth = track_steps;
                    int direction = 1;
                    struct vector *forward = solution->track_plus_motions;
                    if (search_depth < 0) {
                        search_depth = -search_depth;
                        direction = -1;
                        forward = solution->track_minus_motions;
                    }
                    struct vector p = track;
                    for (int i = 0; i < search_depth; ++i) {
                        if (p.u == position.u && p.v == position.v) {
                            track_steps = -i * direction;
                            break;
                        }
                        uint32_t index;
                        if (!lookup_track(solution, p, &index))
                            break;
                        p.u += forward[index].u;
                        p.v += forward[index].v;
                    }
                }
                while (track_steps > 0) {
                    tape[n++] = 't';
                    track_steps--;
                }
                while (track_steps < 0) {
                    tape[n++] = 'g';
                    track_steps++;
                }
                reset_from = n;
                if (n > tape_length)
                    tape_length = n;
                last_end = n;
            } else {
                tape[n] = decode_instruction(inst.instruction);
                last_end = n + 1;
            }
        }
#if 0
        printf("%4u (%4u): %*s", part.arm_number, tape_length, min_tape, "");
        for (uint32_t j = 0; j < tape_length; ++j) {
            if (!tape[j])
                printf(" ");
            else
                printf("%c", tape[j]);
        }
        printf("\n");
#endif
        solution->arm_tape_length[arm_index] = tape_length;
        if (tape_length > solution->tape_period)
            solution->tape_period = tape_length;
        arm_index++;
    }
    // adjust the tape start cycle so everything starts on cycle 1.
    uint32_t tape_start_cycle = UINT32_MAX;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        if (solution->arm_tape_start_cycle[i] < tape_start_cycle)
            tape_start_cycle = solution->arm_tape_start_cycle[i];
    }
    for (uint32_t i = 0; i < solution->number_of_arms; ++i)
        solution->arm_tape_start_cycle[i] -= tape_start_cycle;

    solution->target_number_of_outputs = 6 * pf->output_scale;
    return true;
}
