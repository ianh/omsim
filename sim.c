#include "sim.h"

#include "parse.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double ATOM_RADIUS = 29./82;
static const double SQRT3_2 = 0.8660254037844386;

struct xy_vector {
    double x;
    double y;
};
static struct xy_vector to_xy(struct vector p)
{
    return (struct xy_vector){
        p.v + p.u * 0.5,
        p.u * SQRT3_2,
    };
}
static double xy_len2(struct xy_vector xy)
{
    return xy.x * xy.x + xy.y * xy.y;
}

// hash table functions -- see appendix.
static atom *insert_atom(struct board *board, struct vector query, const char *collision_reason);
static void ensure_capacity(struct board *board, uint32_t potential_insertions);
static void schedule_flag_reset_if_needed(struct board *board, atom *a);

struct vector mechanism_relative_position(struct mechanism m, int32_t du, int32_t dv, int32_t w)
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

static bool conversion_output(struct board *board, bool active, struct mechanism m, int32_t du, int32_t dv)
{
    assert(m.type & CONVERSION_GLYPH);
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    if (board->cycle == 0)
        mark_used_area(board, pos);
    atom *output = lookup_atom(board, pos);
    if ((*output & VALID) && !(*output & REMOVED)) {
        if (board->half_cycle == 2) {
            // conversion glyph outputs appear in the second half-cycle.
            *output &= ~BEING_PRODUCED;
        } else if (active && (*output & BEING_PRODUCED)) {
            report_collision(board, pos, "two conversion glyphs outputting to the same point");
            return true;
        } else if (active && (*output & VAN_BERLO_ATOM)) {
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
    if (board->cycle == 0)
        mark_used_area(board, pos);
    atom *a = lookup_atom(board, pos);
    if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED))
        return (atom *)&empty;
    // only duplication glyphs can see the van berlo's wheel.
    if (!(m.type & (DUPLICATION | VAN_BERLO)) && (*a & VAN_BERLO_ATOM))
        return (atom *)&empty;
    // conversion glyphs can't see bonded or grabbed inputs.
    if ((m.type & CONVERSION_GLYPH) && ((*a & ALL_BONDS) || (*a & GRABBED)))
        return (atom *)&empty;
    return a;
}

static inline void remove_atom(struct board *board, atom *a)
{
    *a |= REMOVED;
}

static inline void transform_atom(atom *a, atom new_type)
{
    *a &= ~ANY_ATOM;
    *a |= new_type;
}

int direction_for_offset(struct vector d)
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

atom bond_direction(struct mechanism m, int32_t du, int32_t dv)
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

// rotate an atom's bonds by rotating the bits which represent those bonds.
static void rotate_bonds(atom *a, int rotation)
{
    rotation = normalize_direction(rotation);
    // first, shift the bond bits by the rotation amount.
    atom bonds = *a & ALL_BONDS;
    bonds <<= rotation;
    // take the overflow bits that were shifted off the end...
    atom mask = 0x3FULL >> (6 - rotation);
    atom overflow = bonds & ((mask * BOND_LOW_BITS) << 6);
    bonds &= ~overflow;
    // ...and shift them back around to the other side.
    bonds |= overflow >> 6;
    *a &= ~ALL_BONDS;
    *a |= bonds;
}

struct vector v_offset_for_direction(int direction)
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

struct vector u_offset_for_direction(int direction)
{
    return v_offset_for_direction(direction + 1);
}

static void apply_conduit(struct solution *solution, struct board *board, struct mechanism m)
{
    struct conduit *conduit = &solution->conduits[m.conduit_index];
    ensure_capacity(board, conduit->number_of_positions);
    struct mechanism other_side = solution->glyphs[conduit->other_side_glyph_index];
    int rotation = direction_for_offset(other_side.direction_v) - direction_for_offset(m.direction_v);
    uint32_t base = 0;
    for (uint32_t j = 0; j < conduit->number_of_molecules; ++j) {
        uint32_t length = conduit->molecule_lengths[j];
        bool valid = true;
        bool consume = true;
        if (board->half_cycle == 1) {
            // if any of the atoms in this molecule have already been consumed
            // by some other glyph, then remove the molecule from the conduit.
            for (uint32_t k = 0; k < length; ++k) {
                struct vector p = conduit->atoms[base + k].position;
                atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
                if (!(*a & VALID) || (*a & REMOVED)) {
                    valid = false;
                    break;
                }
                // mark the shape of the molecule on the board.
                *a |= CONDUIT_SHAPE;
            }
            if (valid) {
                // if a bond has been made between any atom of the molecule and
                // an atom outside the molecule, then don't consume the atoms
                // (but still produce the stored atoms in the second
                // half-cycle).
                for (uint32_t k = 0; k < length && consume; ++k) {
                    struct vector p = conduit->atoms[base + k].position;
                    atom a = *lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
                    for (int bond_direction = 0; bond_direction < 6 && consume; ++bond_direction) {
                        if (!(a & (BOND_LOW_BITS << normalize_direction(bond_direction + direction_for_offset(m.direction_v)))))
                            continue;
                        struct vector d = v_offset_for_direction(bond_direction);
                        atom *b = lookup_atom(board, mechanism_relative_position(m, p.u + d.u, p.v + d.v, 1));
                        if (!(*b & VALID) || (*b & REMOVED) || !(*b & CONDUIT_SHAPE))
                            consume = false;
                    }
                }
            }
            // clear the marked molecule shape.
            for (uint32_t k = 0; k < length; ++k) {
                struct vector p = conduit->atoms[base + k].position;
                atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
                *a &= ~CONDUIT_SHAPE;
            }
            if (!valid) {
                // remove the molecule from the conduit.
                memmove(conduit->atoms + base, conduit->atoms + base + length, (conduit->number_of_positions - base - length) * sizeof(struct atom_at_position));
                memmove(conduit->molecule_lengths + j, conduit->molecule_lengths + j + 1, (conduit->number_of_positions - j - 1) * sizeof(uint32_t));
                j--;
                continue;
            }
        }
        ensure_capacity(board, length);
        for (uint32_t k = 0; k < length; ++k) {
            atom a = conduit->atoms[base + k].atom;
            struct vector delta = conduit->atoms[base + k].position;
            if (board->half_cycle == 1) {
                atom *b = lookup_atom(board, mechanism_relative_position(m, delta.u, delta.v, 1));
                if (consume) {
                    conduit->atoms[base + k].atom = *b;
                    remove_atom(board, b);
                } else {
                    // bonds are consumed even if the atoms aren't.
                    *b &= ~(ALL_BONDS & ~RECENT_BONDS & a);
                }
            } else {
                struct vector p = mechanism_relative_position(other_side, delta.u, delta.v, 1);
                rotate_bonds(&a, rotation);
                *insert_atom(board, p, "conduit output") = a;
            }
        }
        base += length;
    }
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
            bool active = *a & *b & SALT;
            bool c = conversion_output(board, active, m, 1, 0);
            bool d = conversion_output(board, active, m, -1, 1);
            if (c && d && active) {
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
            if (metal && !(*q & ALL_BONDS) && !(*q & GRABBED) && (*q & QUICKSILVER)) {
                remove_atom(board, q);
                transform_atom(a, metal >> 1);
            }
            break;
        }
        case DISPERSION: {
            atom *a = get_atom(board, m, 0, 0);
            bool active = *a & QUINTESSENCE;
            bool b = conversion_output(board, active, m, 0, 1);
            bool c = conversion_output(board, active, m, -1, 1);
            bool d = conversion_output(board, active, m, -1, 0);
            bool e = conversion_output(board, active, m, 0, -1);
            if (b && c && d && e && active) {
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
            atom metal = *a & *b & ANY_METAL & ~GOLD;
            bool c = conversion_output(board, metal, m, 1, 0);
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
            bool active = ((*a | *b | *c | *d) & ANY_ELEMENTAL) == ANY_ELEMENTAL;
            bool e = conversion_output(board, active, m, 0, 0);
            if (e && active) {
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
                atom ab = bond_direction(m, 0, 1);
                atom ba = bond_direction(m, 0, -1);
                // record the bond in the RECENT_BONDS bitfield to prevent the
                // atoms from being consumed during this half-cycle.
                if (*a & ab) {
                    *a &= ~ab;
                    *a |= ab & RECENT_BONDS;
                    schedule_flag_reset_if_needed(board, a);
                }
                if (*b & ba) {
                    *b &= ~ba;
                    *b |= ba & RECENT_BONDS;
                    schedule_flag_reset_if_needed(board, b);
                }
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
            if (board->cycle == 0) {
                mark_used_area(board, mechanism_relative_position(m, 0, 1, 1));
                mark_used_area(board, mechanism_relative_position(m, 1, 0, 1));
                mark_used_area(board, mechanism_relative_position(m, 1, -1, 1));
                mark_used_area(board, mechanism_relative_position(m, 0, -1, 1));
                mark_used_area(board, mechanism_relative_position(m, -1, 0, 1));
                mark_used_area(board, mechanism_relative_position(m, -1, 1, 1));
            }
            if (*a && !(*a & ALL_BONDS) && !(*a & GRABBED))
                remove_atom(board, a);
            break;
        }
        case CONDUIT:
            apply_conduit(solution, board, m);
            break;
        case EQUILIBRIUM:
            if (board->cycle == 0)
                mark_used_area(board, mechanism_relative_position(m, 0, 0, 1));
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

static bool swing_intersection(double swing_radius_squared, struct xy_vector center, struct xy_vector *hit)
{
    double r2 = swing_radius_squared;
    double d2 = xy_len2(center);
    double x = d2 - 4 * ATOM_RADIUS * ATOM_RADIUS + r2;
    if (x * x > 4 * r2 * d2)
        return false;
    double y = sqrt(4 * r2 * d2 - x * x);
    // note that, if there are two intersection points, this chooses the point
    // where the arc *leaves* the circle (the most counterclockwise of the two),
    // which guarantees forward progress.
    hit->x = (x * center.x - y * center.y) / (2 * d2);
    hit->y = (x * center.y + y * center.x) / (2 * d2);
    return true;
}

// returns a positive number if a -> b is a counter-clockwise rotation,
// a negative number if a -> b is a clockwise rotation, and zero if the rotation
// direction is indeterminate (if the points are either equal or antipodal).
static double ccw(struct xy_vector a, struct xy_vector b)
{
    return a.x * b.y - a.y * b.x;
}

static void record_swing_area(struct board *board, struct vector position, struct vector base, int rotation)
{
    struct vector p = position;
    p.u -= base.u;
    p.v -= base.v;
    if (rotation == -1)
        p = (struct vector){ -p.v, p.u + p.v };
    struct xy_vector current = to_xy(p);
    struct xy_vector end = to_xy((struct vector){ p.u + p.v, -p.u });
    double r2 = xy_len2(current);
    if (r2 < 1)
        return;
    while (true) {
        // convert back to grid coordinates.
        int32_t cell_u = (int32_t)floor(current.y / SQRT3_2);
        int32_t cell_v = (int32_t)floor(current.x - 0.5 * current.y / SQRT3_2);
        struct xy_vector min = { 0 };
        struct vector min_cell = { 0 };
        // find the intersection point for the hex at each neighboring grid
        // cell, using the ccw() function to determine which hex the swing arc
        // will encounter next.
        for (int32_t u = cell_u - 1; u <= cell_u + 2; ++u) {
            for (int32_t v = cell_v - 1; v <= cell_v + 2; ++v) {
                struct xy_vector center = to_xy((struct vector){ u, v });
                struct xy_vector hit;
                if (!swing_intersection(r2, center, &hit))
                    continue;
                if (ccw(current, hit) <= 0)
                    continue;
                if ((min.x != 0 || min.y != 0) && ccw(hit, min) <= 0)
                    continue;
                min = hit;
                min_cell.u = u;
                min_cell.v = v;
            }
        }
        assert(min.x != 0 || min.y != 0);
        if (ccw(min, end) <= 0)
            break;
        mark_used_area(board, (struct vector){
            base.u + min_cell.u,
            base.v + min_cell.v,
        });
        // advance the current point to the point of encounter.
        current = min;
    }
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
            if (!(m.atom & (BOND_LOW_BITS << bond_direction) & ~RECENT_BONDS))
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
    if (solution->tape_period == 0)
        return;
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
            fprintf(stderr, "trying to extend/retract a non-piston arm");
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
                fprintf(stderr, "trying to move an arm along a track that isn't on a track at %" PRId32 " %" PRId32 "\n", m->position.u, m->position.v);
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
            struct vector saved_u = m->direction_u;
            struct vector saved_v = m->direction_v;
            while (true) {
                // xx do area somewhere else?
                struct vector p = mechanism_relative_position(*m, offset.u, offset.v, 1);
                mark_used_area(board, p);
                if (inst == 'a')
                    record_swing_area(board, p, m->position, 1);
                if (inst == 'd')
                    record_swing_area(board, p, m->position, -1);
                if (vectors_equal(m->direction_u, zero_vector))
                    break;
                adjust_axis_magnitude(&m->direction_u, -1);
                adjust_axis_magnitude(&m->direction_v, -1);
            }
            m->direction_u = saved_u;
            m->direction_v = saved_v;
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
                // allow conduits to transport this atom (and any atoms bonded to it).
                *a |= BEING_DROPPED;
                schedule_flag_reset_if_needed(board, a);
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
        if (m.rotation) {
            rotate_bonds(&m.atom, m.rotation);
            record_swing_area(board, m.position, m.base, m.rotation);
        }
        atom *a = insert_atom(board, to, "atom moved on top of another atom");
        *a = VALID | m.atom;

        // xx collision detection / error handling
    }
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        // xx definitely do this in a cleaner way...
        struct mechanism *m = &solution->arms[i];
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
            struct vector saved_u = m->direction_u;
            struct vector saved_v = m->direction_v;
            while (true) {
                struct vector p = mechanism_relative_position(*m, offset.u, offset.v, 1);
                mark_used_area(board, p);
                if (vectors_equal(m->direction_u, zero_vector))
                    break;
                adjust_axis_magnitude(&m->direction_u, -1);
                adjust_axis_magnitude(&m->direction_v, -1);
            }
            m->direction_u = saved_u;
            m->direction_v = saved_v;
        }
    }
}

static bool fill_conduit_molecule(struct solution *solution, struct board *board, struct conduit *conduit, uint32_t *atom_next, uint32_t *atom_cursor)
{
    struct mechanism m = solution->glyphs[conduit->glyph_index];
    int rotation = direction_for_offset(m.direction_v);
    while (*atom_cursor < *atom_next) {
        struct atom_at_position ap = conduit->atoms[*atom_cursor];
        for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
            if (!(ap.atom & (BOND_LOW_BITS << normalize_direction(bond_direction + rotation)) & ~RECENT_BONDS))
                continue;
            struct vector p = ap.position;
            struct vector d = v_offset_for_direction(bond_direction);
            p.u += d.u;
            p.v += d.v;
            atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
            if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED) || (*a & VISITED))
                continue;
            if (!(*a & CONDUIT_SHAPE) || (*a & GRABBED))
                return false;
            conduit->atoms[*atom_next].atom = *a & ~BEING_DROPPED & ~CONDUIT_SHAPE;
            conduit->atoms[*atom_next].position = p;
            (*atom_next)++;
            *a |= VISITED;
        }
        (*atom_cursor)++;
    }
    return true;
}

static void fill_conduits(struct solution *solution, struct board *board)
{
    if (board->half_cycle != 1)
        return;
    for (uint32_t i = 0; i < solution->number_of_conduits; ++i) {
        struct conduit *conduit = &solution->conduits[i];
        struct mechanism m = solution->glyphs[conduit->glyph_index];
        // mark every atom on the conduit with the CONDUIT_SHAPE flag.
        for (uint32_t j = 0; j < conduit->number_of_positions; ++j) {
            struct vector p = conduit->positions[j];
            atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
            if ((*a & VALID) && !(*a & REMOVED))
                *a |= CONDUIT_SHAPE;
        }
        // discover the molecules attached to any atoms dropped on the conduit.
        // if the molecule is completely inside the conduit, add its atoms to
        // conduit->atoms and record its length in conduit->molecule_lengths.
        conduit->number_of_molecules = 0;
        uint32_t atom_next = 0;
        uint32_t atom_cursor = 0;
        for (uint32_t j = 0; j < conduit->number_of_positions; ++j) {
            struct vector p = conduit->positions[j];
            atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
            if ((*a & VALID) && !(*a & REMOVED) && !(*a & BEING_PRODUCED) && !(*a & VISITED) && (*a & BEING_DROPPED) && !(*a & GRABBED)) {
                conduit->atoms[atom_next].atom = *a & ~BEING_DROPPED & ~CONDUIT_SHAPE;
                conduit->atoms[atom_next].position = p;
                *a |= VISITED;
                atom_next++;
                uint32_t saved_cursor = atom_cursor;
                if (fill_conduit_molecule(solution, board, conduit, &atom_next, &atom_cursor)) {
                    // the molecule matches the conduit shape; clear the
                    // BEING_DROPPED flag so no other conduit can see it.
                    for (uint32_t k = saved_cursor; k < atom_cursor; ++k) {
                        struct vector p = conduit->atoms[k].position;
                        p = mechanism_relative_position(m, p.u, p.v, 1);
                        *lookup_atom(board, p) &= ~BEING_DROPPED;
                    }
                    conduit->molecule_lengths[conduit->number_of_molecules++] = atom_cursor - saved_cursor;
                } else {
                    // the molecule doesn't match the conduit shape.
                    atom_next = saved_cursor;
                    atom_cursor = saved_cursor;
                }
            }
        }
        // clear the CONDUIT_SHAPE and VISITED flags so the next conduit can use
        // them.
        for (uint32_t j = 0; j < conduit->number_of_positions; ++j) {
            struct vector p = conduit->positions[j];
            atom *a = lookup_atom(board, mechanism_relative_position(m, p.u, p.v, 1));
            *a &= ~CONDUIT_SHAPE & ~VISITED;
        }
    }
}

static void flag_blocked_inputs(struct solution *solution, struct board *board)
{
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i)
        solution->inputs_and_outputs[i].type &= ~BLOCKED;
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        if (!(io->type & INPUT))
            continue;
        for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
            atom a = *lookup_atom(board, io->atoms[j].position);
            if ((a & VALID) && !(a & REMOVED) && !(a & BEING_PRODUCED)) {
                io->type |= BLOCKED;
                break;
            }
        }
    }
}

static void spawn_inputs(struct solution *solution, struct board *board)
{
    size_t active = board->active_input_or_output;
    size_t i = active == UINT32_MAX ? 0 : active;
    for (; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        if (!(io->type & INPUT) || (io->type & BLOCKED))
            continue;
        if (i != active && (io->type & INTERRUPT)) {
            board->active_input_or_output = i;
            return;
        }
        uint32_t n = io->number_of_atoms;
        ensure_capacity(board, n);
        for (uint32_t j = 0; j < n; ++j) {
            atom input = io->atoms[j].atom;
            struct vector p = io->atoms[j].position;
            *insert_atom(board, p, "overlapped inputs") = VALID | input;
        }
    }
    board->active_input_or_output = UINT32_MAX;
}

static void consume_outputs(struct solution *solution, struct board *board)
{
    size_t active = board->active_input_or_output;
    size_t i = active == UINT32_MAX ? 0 : active;
    for (; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        if (!(io->type & OUTPUT))
            continue;
        // first, check the entire output to see if it matches.
        bool match = true;
        for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
            atom output = io->atoms[j].atom;
            if (output & REPEATING_OUTPUT_PLACEHOLDER)
                continue;
            atom *a = lookup_atom(board, io->atoms[j].position);
            // printf("checking %d %d... ", io->atoms[j].position.u, io->atoms[j].position.v);
            if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED)) {
                // printf("no atom\n");
                match = false;
                break;
            }
            // variable outputs match any atom.
            if (output & VARIABLE_OUTPUT) {
                output &= ~VARIABLE_OUTPUT;
                output |= *a & ANY_ATOM;
            }
            if (io->type & REPEATING_OUTPUT) {
                if ((*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS & output))) != output) {
                    // printf("did not match: %llx vs %llx\n", (*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS & output))), output);
                    match = false;
                    break;
                }
            } else {
                if ((*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS))) != output || (*a & GRABBED)) {
                    // printf("did not match: %llx vs %llx\n", (*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS))), output);
                    match = false;
                    break;
                }
            }
            // printf("match!\n");
        }
        // if the output is a match, first trigger an interrupt if necessary.
        // then remove the output and increment the output counter.
        if (match) {
            if (i != active && (io->type & INTERRUPT)) {
                board->active_input_or_output = i;
                return;
            }
            if (io->type & REPEATING_OUTPUT)
                io->number_of_outputs = REPEATING_OUTPUT_REPETITIONS;
            else {
                for (uint32_t j = 0; j < io->number_of_atoms; ++j)
                    remove_atom(board, lookup_atom(board, io->atoms[j].position));
                io->number_of_outputs++;
            }
        }
    }
    board->active_input_or_output = UINT32_MAX;
}

static void reset_temporary_flags(struct board *board)
{
    for (uint32_t i = 0; i < board->flag_reset_length; ++i)
        *board->flag_reset[i] &= ~TEMPORARY_FLAGS;
    board->flag_reset_length = 0;
}

static bool check_completion(struct solution *solution)
{
    uint64_t min = UINT64_MAX;
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        if (!(solution->inputs_and_outputs[i].type & OUTPUT))
            continue;
        if ((solution->inputs_and_outputs[i].type & REPEATING_OUTPUT) &&
         solution->inputs_and_outputs[i].number_of_outputs == REPEATING_OUTPUT_REPETITIONS)
            continue;
        uint64_t count = solution->inputs_and_outputs[i].number_of_outputs;
        if (count < min)
            min = count;
    }
    return min >= solution->target_number_of_outputs;
}

enum run_result run(struct solution *solution, struct board *board)
{
    for (; board->half_cycle <= 2; board->half_cycle++) {
        uint32_t io = board->active_input_or_output;
        if (io != UINT32_MAX) {
            if (solution->inputs_and_outputs[io].type & INPUT)
                goto continue_with_inputs;
            else
                goto continue_with_outputs;
        }
        perform_arm_instructions(solution, board);
        fill_conduits(solution, board);
        flag_blocked_inputs(solution, board);
continue_with_inputs:
        spawn_inputs(solution, board);
        if (board->active_input_or_output != UINT32_MAX)
            return INPUT_OUTPUT;
        reset_temporary_flags(board);
        apply_glyphs(solution, board);
continue_with_outputs:
        consume_outputs(solution, board);
        if (board->active_input_or_output != UINT32_MAX)
            return INPUT_OUTPUT;
    }
    board->complete = check_completion(solution);
    board->cycle++;
    board->half_cycle = 1;
    return FINISHED_CYCLE;
}

static void create_van_berlo_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv, atom element)
{
    ensure_capacity(board, 1);
    struct vector p = mechanism_relative_position(m, du, dv, 1);
    atom *a = insert_atom(board, p, "van berlo overlap");
    *a = VALID | GRABBED_ONCE | VAN_BERLO_ATOM | element;
}

void initial_setup(struct solution *solution, struct board *board, uint32_t intial_board_size)
{
    board->half_cycle = 1;
    board->active_input_or_output = UINT32_MAX;
    // in order for lookups to work, the hash table has to be allocated using
    // ensure_capacity().
    if (intial_board_size < 1)
        intial_board_size = 1;
    ensure_capacity(board, intial_board_size);
    for (uint32_t i = 0; i < solution->track_table_size; ++i) {
        // xx do this in a cleaner way?
        struct vector p = solution->track_positions[i];
        if (p.u == INT32_MIN && p.v == INT32_MIN)
            continue;
        mark_used_area(board, p);
    }
    for (uint32_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        for (uint32_t j = 0; j < io->number_of_atoms; ++j)
            mark_used_area(board, io->atoms[j].position);
    }

    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        if (!(solution->arms[i].type & VAN_BERLO))
            continue;
        solution->arms[i].type |= GRABBING;
        solution->arms[i].type |= GRABBING_EVERYTHING;
        create_van_berlo_atom(board, solution->arms[i], 0, 1, SALT);
        create_van_berlo_atom(board, solution->arms[i], 1, 0, WATER);
        create_van_berlo_atom(board, solution->arms[i], 1, -1, AIR);
        create_van_berlo_atom(board, solution->arms[i], 0, -1, SALT);
        create_van_berlo_atom(board, solution->arms[i], -1, 0, FIRE);
        create_van_berlo_atom(board, solution->arms[i], -1, 1, EARTH);
    }
    // van berlo's wheel can block inputs.
    flag_blocked_inputs(solution, board);
    spawn_inputs(solution, board);
}

void destroy(struct solution *solution, struct board *board)
{
    free(solution->glyphs);
    free(solution->arms);
    for (uint32_t i = 0; i < solution->number_of_arms; ++i)
        free(solution->arm_tape[i]);
    free(solution->arm_tape);
    free(solution->arm_tape_length);
    free(solution->arm_tape_start_cycle);
    free(solution->track_positions);
    free(solution->track_plus_motions);
    free(solution->track_minus_motions);
    for (size_t i = 0; i < solution->number_of_conduits; ++i) {
        free(solution->conduits[i].positions);
        free(solution->conduits[i].atoms);
        free(solution->conduits[i].molecule_lengths);
    }
    free(solution->conduits);
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i)
        free(solution->inputs_and_outputs[i].atoms);
    free(solution->inputs_and_outputs);

    free(board->atoms_at_positions);
    free(board->flag_reset);
    free(board->movements.elements);
}

uint32_t used_area(struct board *board)
{
    // xx record/replay swings for extra efficiency
    return board->used;
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

bool lookup_track(struct solution *solution, struct vector query, uint32_t *index)
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
        if (vectors_equal(p, query))
            return true;
        *index = (*index + 1) & mask;
        if (*index == (hash & mask))
            abort();
    }
}

static struct atom_at_position *lookup_atom_at_position(struct board *board, struct vector query)
{
    uint32_t hash = fnv(&query, sizeof(query));
    uint32_t mask = board->capacity - 1;
    uint32_t index = hash & mask;
    while (true) {
        struct atom_at_position *a = &board->atoms_at_positions[index];
        if (!(a->atom & VALID))
            return a;
        if (vectors_equal(a->position, query))
            return a;
        index = (index + 1) & mask;
        if (index == (hash & mask))
            abort();
    }
}

atom *lookup_atom(struct board *board, struct vector query)
{
    return &lookup_atom_at_position(board, query)->atom;
}

void mark_used_area(struct board *board, struct vector point)
{
    ensure_capacity(board, 1);
    struct atom_at_position *a = lookup_atom_at_position(board, point);
    if (a->atom & VALID)
        return;
    a->position = point;
    a->atom = VALID | REMOVED;
    board->used++;
}

static atom *insert_atom(struct board *board, struct vector query, const char *collision_reason)
{
    struct atom_at_position *a = lookup_atom_at_position(board, query);
    if ((a->atom & VALID) && !(a->atom & REMOVED))
        report_collision(board, query, collision_reason);
    a->position = query;
    if (!(a->atom & REMOVED))
        board->used++;
    return &a->atom;
}

static void ensure_capacity(struct board *board, uint32_t potential_insertions)
{
    uint32_t n = board->capacity;
    if (n == 0)
        n = 16;
    while (2 * n <= 7 * (board->used + potential_insertions))
        n *= 2;
    if (n == board->capacity)
        return;
    struct board old = *board;
    board->atoms_at_positions = calloc(n, sizeof(*board->atoms_at_positions));
    board->capacity = n;
    board->used = 0;
    board->flag_reset_length = 0;
    for (uint32_t i = 0; i < old.capacity; ++i) {
        atom a = old.atoms_at_positions[i].atom;
        if (!(a & VALID))
            continue;
        struct vector position = old.atoms_at_positions[i].position;
        atom *b = insert_atom(board, position, "reinserting colliding atoms");
        *b = a;
        schedule_flag_reset_if_needed(board, b);
    }
    free(old.atoms_at_positions);
}

static void schedule_flag_reset_if_needed(struct board *board, atom *a)
{
    if (!(*a & TEMPORARY_FLAGS))
        return;
    if (board->flag_reset_length >= board->flag_reset_capacity) {
        uint32_t n = (board->flag_reset_length + 1) * 2;
        atom **flag_reset = realloc(board->flag_reset, n * sizeof(atom *));
        if (!flag_reset)
            abort();
        board->flag_reset = flag_reset;
        board->flag_reset_capacity = n;
    }
    board->flag_reset[board->flag_reset_length++] = a;
}
