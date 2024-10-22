#include "sim.h"

#include "collision.h"
#include "parse.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct atom_ref_at_position {
    atom *atom;
    struct vector position;
};

// hash table functions -- see appendix.
static void rehash(struct atom_grid *grid, uint32_t size);
static void schedule_flag_reset_if_needed(struct board *board, atom *a);

static atom mark_used_area_with_overlap(struct board *board, struct vector point, uint64_t *overlap);

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
        return false;
    }
    return true;
}

static void produce_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv, atom a)
{
    assert(m.type & CONVERSION_GLYPH);
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    insert_atom(board, pos, VALID | BEING_PRODUCED | a, "conversion glyph output");
}

static struct atom_ref_at_position get_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv)
{
    static const atom empty;
    struct vector pos = mechanism_relative_position(m, du, dv, 1);
    // conversion glyphs don't consume any inputs in the second half-cycle.
    if ((m.type & CONVERSION_GLYPH) && board->half_cycle == 2)
        return (struct atom_ref_at_position){ (atom *)&empty, pos };
    atom *a = lookup_atom(board, pos);
    if (!(*a & VALID))
        return (struct atom_ref_at_position){ (atom *)&empty, pos };
    // arms need to see removed-but-moving atoms to perform simultaneous motion checks.
    if ((m.type & ANY_ARM) && (*a & REMOVED) && (*a & MOVED))
        return (struct atom_ref_at_position){ a, pos };
    if ((*a & REMOVED) || ((*a & BEING_PRODUCED) && !(m.type & UNBONDING)))
        return (struct atom_ref_at_position){ (atom *)&empty, pos };
    // only duplication glyphs can see the van berlo's wheel.
    if (!(m.type & (DUPLICATION | VAN_BERLO)) && (*a & VAN_BERLO_ATOM))
        return (struct atom_ref_at_position){ (atom *)&empty, pos };
    // conversion glyphs can't see bonded or grabbed inputs.
    if ((m.type & CONVERSION_GLYPH) && ((*a & ALL_BONDS) || (*a & GRABBED)))
        return (struct atom_ref_at_position){ (atom *)&empty, pos };
    return (struct atom_ref_at_position){ a, pos };
}

__attribute__((noinline))
static void remove_overlapping_atom(struct board *board, struct atom_ref_at_position a)
{
    for (uint32_t i = 0; i < board->number_of_overlapped_atoms; ++i) {
        if (vectors_equal(board->overlapped_atoms[i].position, a.position)) {
            *a.atom = board->overlapped_atoms[i].atom;
            schedule_flag_reset_if_needed(board, a.atom);
            memmove(board->overlapped_atoms + i,
                    board->overlapped_atoms + i + 1,
                    (board->number_of_overlapped_atoms - i - 1) * sizeof(struct atom_at_position));
            board->number_of_overlapped_atoms--;
            return;
        }
    }
}

static inline bool remove_atom(struct board *board, struct atom_ref_at_position a)
{
    if (*a.atom & OVERLAPS_ATOMS) {
        remove_overlapping_atom(board, a);
        return false;
    } else {
        *a.atom = VALID | REMOVED;
        return true;
    }
}

static inline void transform_atom(atom *a, atom new_type)
{
    *a &= ~ANY_ATOM;
    *a |= new_type;
}

int direction_for_offset(struct vector d)
{
    if (d.u == 1 && d.v == 0)
        return 0;
    else if (d.u == 0 && d.v == 1)
        return 1;
    else if (d.u == -1 && d.v == 1)
        return 2;
    else if (d.u == -1 && d.v == 0)
        return 3;
    else if (d.u == 0 && d.v == -1)
        return 4;
    else if (d.u == 1 && d.v == -1)
        return 5;
    return -1;
}

int angular_distance_between_grabbers(enum mechanism_type type)
{
    switch (type & ANY_ARM) {
    case ARM:
    case PISTON:
        return 6;
    case TWO_ARM:
        return 3;
    case THREE_ARM:
        return 2;
    case SIX_ARM:
    case VAN_BERLO:
    default:
        return 1;
    }
}

atom bond_direction(struct mechanism m, int32_t du, int32_t dv)
{
    atom base = BOND_LOW_BITS;
    int dir = direction_for_offset(mechanism_relative_position(m, du, dv, 0));
    if (dir < 0) {
        // fprintf(stderr, "internal error: non-orthonormal bond direction (scaled bonder?)\n");
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
    if (rotation == 0)
        return;
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

struct vector u_offset_for_direction(int direction)
{
    switch (normalize_direction(direction)) {
    case 0: return (struct vector){1, 0};
    case 1: return (struct vector){0, 1};
    case 2: return (struct vector){-1, 1};
    case 3: return (struct vector){-1, 0};
    case 4: return (struct vector){0, -1};
    case 5: return (struct vector){1, -1};
    default: abort();
    }
}

struct vector v_offset_for_direction(int direction)
{
    return u_offset_for_direction(direction + 1);
}

static void apply_conduit(struct solution *solution, struct board *board, struct mechanism m)
{
    struct conduit *conduit = &solution->conduits[m.conduit_index];
    struct mechanism other_side = solution->glyphs[conduit->other_side_glyph_index];
    uint32_t base = 0;

    if (board->half_cycle == 1) {
        int rotation = direction_for_offset(other_side.direction_u) - direction_for_offset(m.direction_u);
        for (uint32_t j = 0; j < conduit->number_of_molecules; ++j) {
            uint32_t length = conduit->molecule_lengths[j];
            bool valid = true;
            bool consume = true;
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
                        if (!(a & (BOND_LOW_BITS << normalize_direction(bond_direction + direction_for_offset(m.direction_u)))))
                            continue;
                        struct vector d = u_offset_for_direction(bond_direction);
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
                memmove(conduit->molecule_lengths + j, conduit->molecule_lengths + j + 1, (conduit->number_of_molecules - j - 1) * sizeof(uint32_t));
                conduit->number_of_molecules--;
                j--;
                continue;
            }
            for (uint32_t k = 0; k < length; ++k) {
                struct vector delta = conduit->atoms[base + k].position;
                struct vector p = mechanism_relative_position(m, delta.u, delta.v, 1);
                atom a = conduit->atoms[base + k].atom;
                atom *b = lookup_atom(board, p);
                if (consume) {
                    conduit->atoms[base + k].atom = *b;
                    remove_atom(board, (struct atom_ref_at_position){ b, p });
                } else {
                    // bonds are consumed even if the atoms aren't.
                    *b &= ~(ALL_BONDS & ~RECENT_BONDS & a);
                }
                p = mechanism_relative_position(other_side, delta.u, delta.v, 1);
                rotate_bonds(&a, rotation);
                insert_atom(board, p, a | BEING_PRODUCED, "conduit output");
            }
            base += length;
        }
    } else {
        conduit = &solution->conduits[other_side.conduit_index];
        for (uint32_t j = 0; j < conduit->number_of_molecules; ++j) {
            uint32_t length = conduit->molecule_lengths[j];
            for (uint32_t k = 0; k < length; ++k) {
                struct vector delta = conduit->atoms[base + k].position;
                struct vector p = mechanism_relative_position(m, delta.u, delta.v, 1);
                *lookup_atom(board, p) &= ~BEING_PRODUCED;
            }
            base += length;
        }
    }
}

__attribute__((noinline))
static void unbond_overlapping_atoms(struct board *board, struct atom_ref_at_position a, struct atom_ref_at_position b, atom ab, atom ba)
{
    // this isn't actually correct, but just remove all relevant bonds from all overlapping atoms.
    for (uint32_t i = 0; i < board->number_of_overlapped_atoms; ++i) {
        struct atom_at_position overlap = board->overlapped_atoms[i];
        if (vectors_equal(overlap.position, a.position) && (overlap.atom & ab)) {
            board->overlapped_atoms[i].atom &= ~ab;
            board->overlapped_atoms[i].atom |= ab & RECENT_BONDS;
        }
        if (vectors_equal(overlap.position, b.position) && (overlap.atom & ba)) {
            board->overlapped_atoms[i].atom &= ~ba;
            board->overlapped_atoms[i].atom |= ba & RECENT_BONDS;
        }
    }
}

static void apply_glyphs(struct solution *solution, struct board *board)
{
    size_t n = solution->number_of_glyphs;
    for (size_t i = 0; i < n; ++i) {
        struct mechanism m = solution->glyphs[i];
        switch (m.type & ANY_GLYPH) {
        case CALCIFICATION: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            if (*a.atom & ANY_ELEMENTAL)
                transform_atom(a.atom, SALT);
            break;
        }
        case ANIMISMUS: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            struct atom_ref_at_position b = get_atom(board, m, 1, 0);
            bool active = *a.atom & *b.atom & SALT;
            bool c = conversion_output(board, active, m, 0, 1);
            bool d = conversion_output(board, active, m, 1, -1);
            if (c && d && active) {
                remove_atom(board, a);
                remove_atom(board, b);
                produce_atom(board, m, 0, 1, VITAE);
                produce_atom(board, m, 1, -1, MORS);
            }
            break;
        }
        case PROJECTION: {
            struct atom_ref_at_position q = get_atom(board, m, 0, 0);
            struct atom_ref_at_position a = get_atom(board, m, 1, 0);
            atom metal = *a.atom & ANY_METAL & ~GOLD;
            if (metal && !(*q.atom & ALL_BONDS) && !(*q.atom & GRABBED) && (*q.atom & QUICKSILVER)) {
                remove_atom(board, q);
                transform_atom(a.atom, metal >> 1);
            }
            break;
        }
        case DISPERSION: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            bool active = *a.atom & QUINTESSENCE;
            bool b = conversion_output(board, active, m, 1, 0);
            bool c = conversion_output(board, active, m, 1, -1);
            bool d = conversion_output(board, active, m, 0, -1);
            bool e = conversion_output(board, active, m, -1, 0);
            if (b && c && d && e && active) {
                remove_atom(board, a);
                produce_atom(board, m, 1, 0, EARTH);
                produce_atom(board, m, 1, -1, WATER);
                produce_atom(board, m, 0, -1, FIRE);
                produce_atom(board, m, -1, 0, AIR);
            }
            break;
        }
        case PURIFICATION: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            struct atom_ref_at_position b = get_atom(board, m, 1, 0);
            atom metal = *a.atom & *b.atom & ANY_METAL & ~GOLD;
            bool c = conversion_output(board, metal, m, 0, 1);
            if (c && metal) {
                remove_atom(board, a);
                remove_atom(board, b);
                produce_atom(board, m, 0, 1, metal >> 1);
            }
            break;
        }
        case DUPLICATION: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            struct atom_ref_at_position b = get_atom(board, m, 1, 0);
            atom elemental = *a.atom & ANY_ELEMENTAL;
            if (elemental && (*b.atom & SALT) && !(*b.atom & VAN_BERLO_ATOM))
                transform_atom(b.atom, elemental);
            break;
        }
        case UNIFICATION: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 1);
            struct atom_ref_at_position b = get_atom(board, m, -1, 1);
            struct atom_ref_at_position c = get_atom(board, m, 0, -1);
            struct atom_ref_at_position d = get_atom(board, m, 1, -1);
            bool active = ((*a.atom | *b.atom | *c.atom | *d.atom) & ANY_ELEMENTAL) == ANY_ELEMENTAL;
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
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            struct atom_ref_at_position b = get_atom(board, m, 1, 0);
            if (*a.atom && *b.atom)
                add_bond(m, a.atom, b.atom, 1, 0, NORMAL_BONDS, TRIPLEX_BONDS);
            break;
        }
        case UNBONDING: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            struct atom_ref_at_position b = get_atom(board, m, 1, 0);
            if (*a.atom && *b.atom) {
                atom ab = bond_direction(m, 1, 0);
                atom ba = bond_direction(m, -1, 0);
                // record the bond in the RECENT_BONDS bitfield to prevent the
                // atoms from being consumed during this half-cycle.
                if (*a.atom & ab) {
                    *a.atom &= ~ab;
                    *a.atom |= ab & RECENT_BONDS;
                    schedule_flag_reset_if_needed(board, a.atom);
                }
                if (*b.atom & ba) {
                    *b.atom &= ~ba;
                    *b.atom |= ba & RECENT_BONDS;
                    schedule_flag_reset_if_needed(board, b.atom);
                }
                if ((*a.atom & OVERLAPS_ATOMS) || (*b.atom & OVERLAPS_ATOMS))
                    unbond_overlapping_atoms(board, a, b, ab, ba);
            }
            break;
        }
        case TRIPLEX_BONDING: {
            struct atom_ref_at_position ky = get_atom(board, m, 0, 0);
            struct atom_ref_at_position yr = get_atom(board, m, 0, 1);
            struct atom_ref_at_position rk = get_atom(board, m, 1, 0);
            if (*ky.atom && *yr.atom && (*ky.atom & *yr.atom & FIRE))
                add_bond(m, ky.atom, yr.atom, 0, 1, TRIPLEX_Y_BONDS, NORMAL_BONDS);
            if (*yr.atom && *rk.atom && (*yr.atom & *rk.atom & FIRE))
                add_bond(m, yr.atom, rk.atom, 1, -1, TRIPLEX_R_BONDS, NORMAL_BONDS);
            if (*rk.atom && *ky.atom && (*rk.atom & *ky.atom & FIRE))
                add_bond(m, rk.atom, ky.atom, -1, 0, TRIPLEX_K_BONDS, NORMAL_BONDS);
            break;
        }
        case MULTI_BONDING: {
            struct atom_ref_at_position center = get_atom(board, m, 0, 0);
            struct atom_ref_at_position a = get_atom(board, m, 1, 0);
            struct atom_ref_at_position b = get_atom(board, m, 0, -1);
            struct atom_ref_at_position c = get_atom(board, m, -1, 1);
            if (*center.atom && *a.atom)
                add_bond(m, a.atom, center.atom, -1, 0, NORMAL_BONDS, TRIPLEX_BONDS);
            if (*center.atom && *b.atom)
                add_bond(m, b.atom, center.atom, 0, 1, NORMAL_BONDS, TRIPLEX_BONDS);
            if (*center.atom && *c.atom)
                add_bond(m, c.atom, center.atom, 1, -1, NORMAL_BONDS, TRIPLEX_BONDS);
            break;
        }
        case DISPOSAL: {
            struct atom_ref_at_position a = get_atom(board, m, 0, 0);
            if (*a.atom && !(*a.atom & ALL_BONDS) && !(*a.atom & GRABBED))
                remove_atom(board, a);
            break;
        }
        case CONDUIT:
            apply_conduit(solution, board, m);
            break;
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

static void record_swing_area(struct board *board, struct vector base, struct vector position, int rotation)
{
    struct vector offset = { position.u - base.u, position.v - base.v };
    if (rotation == -1)
        offset = (struct vector){ offset.u + offset.v, -offset.u };
    int32_t arm_length = abs(offset.u) | abs(offset.v);
    offset = normalize_axis(offset);
    switch (arm_length) {
    case 3:
        mark_used_area(board, (struct vector){ base.u + 3 * offset.u - 1 * offset.v, base.v + 1 * offset.u + 4 * offset.v });
        mark_used_area(board, (struct vector){ base.u + 2 * offset.u - 2 * offset.v, base.v + 2 * offset.u + 4 * offset.v });
        mark_used_area(board, (struct vector){ base.u + 1 * offset.u - 3 * offset.v, base.v + 3 * offset.u + 4 * offset.v });
        // fall through
    case 2:
        mark_used_area(board, (struct vector){ base.u + 2 * offset.u - 1 * offset.v, base.v + 1 * offset.u + 3 * offset.v });
        mark_used_area(board, (struct vector){ base.u + 1 * offset.u - 1 * offset.v, base.v + 1 * offset.u + 2 * offset.v });
        mark_used_area(board, (struct vector){ base.u + 1 * offset.u - 2 * offset.v, base.v + 2 * offset.u + 3 * offset.v });
        // fall through
    case 1:
        return;
    default:
        abort();
    }
}

static size_t reserve_movement_index(struct board *board)
{
    if (board->movements.length >= MAX_MOVEMENTS)
        abort();
    size_t capacity = board->movements.capacity;
    while (board->movements.length >= capacity)
        capacity = 4 * (capacity + 12) / 3;
    if (capacity != board->movements.capacity) {
        struct movement *elems = realloc(board->movements.movements,
         sizeof(struct movement) * capacity);
        if (!elems)
            abort();
        board->movements.movements = elems;
        board->movements.capacity = capacity;
    }
    return board->movements.length++;
}

static bool movements_equal(struct movement a, struct movement b)
{
    if (a.rotation == 0 && b.rotation == 0)
        return vectors_equal(a.translation, b.translation);
    if (a.rotation != b.rotation)
        return false;
    struct vector a_center = a.type == SWING_MOVEMENT ? a.base : a.absolute_grab_position;
    struct vector b_center = b.type == SWING_MOVEMENT ? b.base : b.absolute_grab_position;
    return vectors_equal(a_center, b_center);
}

static void move_atom(struct board *board, atom *a, struct vector position, size_t movement_index, uint32_t *chain_atom_list)
{
    size_t capacity = board->moving_atoms.capacity;
    while (board->moving_atoms.length >= capacity)
        capacity = 4 * (capacity + 12) / 3;
    if (capacity != board->moving_atoms.capacity) {
        struct atom_at_position *elems = realloc(board->moving_atoms.atoms_at_positions,
         sizeof(struct atom_at_position) * capacity);
        if (!elems)
            abort();
        board->moving_atoms.atoms_at_positions = elems;
        board->moving_atoms.capacity = capacity;
    }
    board->moving_atoms.atoms_at_positions[board->moving_atoms.length++] = (struct atom_at_position){
        .atom = *a,
        .position = position,
    };
    if (*a & IS_CHAIN_ATOM)
        move_chain_atom_to_list(board, lookup_chain_atom(board, position), chain_atom_list);
    if (remove_atom(board, (struct atom_ref_at_position){ a, position }))
        *a |= MOVED | MOVEMENT_INDEX(movement_index);
}

static void move_atoms(struct board *board, atom *a, struct movement movement)
{
    if ((*a & VALID) && (*a & REMOVED) && (*a & MOVED)) {
        if (!movements_equal(movement, board->movements.movements[GET_MOVEMENT_INDEX(*a)]))
            report_collision(board, movement.absolute_grab_position, "atom moved in two directions simultaneously");
        return;
    }
    movement.first_chain_atom = UINT32_MAX;
    movement.first_atom_index = board->moving_atoms.cursor;
    size_t movement_index = reserve_movement_index(board);
    move_atom(board, a, movement.absolute_grab_position, movement_index, &movement.first_chain_atom);
    // do a breadth-first search over the molecule, removing each atom from the
    // board as it's discovered.
    while (board->moving_atoms.cursor < board->moving_atoms.length) {
        struct atom_at_position m = board->moving_atoms.atoms_at_positions[board->moving_atoms.cursor];
        for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
            if (!(m.atom & (BOND_LOW_BITS << bond_direction) & ~RECENT_BONDS))
                continue;
            struct vector p = m.position;
            struct vector d = u_offset_for_direction(bond_direction);
            p.u += d.u;
            p.v += d.v;
            atom *b = lookup_atom(board, p);
            if (!(*b & VALID) || (*b & REMOVED) || (*b & BEING_PRODUCED))
                continue;
            move_atom(board, b, p, movement_index, &movement.first_chain_atom);
        }
        board->moving_atoms.cursor++;
    }
    movement.number_of_atoms = board->moving_atoms.cursor - movement.first_atom_index;
    board->movements.movements[movement_index] = movement;
}

static void apply_movement_to_position(struct vector base, struct vector translation, struct vector u, struct vector v, struct vector *position)
{
    struct vector delta = *position;
    delta.u -= base.u;
    delta.v -= base.v;
    *position = (struct vector){
        u.u * delta.u + v.u * delta.v + base.u + translation.u,
        u.v * delta.u + v.v * delta.v + base.v + translation.v,
    };
}

static bool range_of_periods_where_motion_is_in_interval(int32_t motion, int32_t interval_min, int32_t interval_max, int32_t *entry_period, int32_t *exit_period)
{
    if (motion == 0) {
        // if there's no motion, either we stay in the interval forever (if zero
        // is in the interval) or we never enter it.
        *entry_period = INT32_MIN;
        *exit_period = INT32_MAX;
        return 0 >= interval_min && 0 <= interval_max;
    } else if (motion < 0) {
        // if the motion is negative, flip everything so it can be positive.
        // note that max and min change roles when their signs change.
        motion = -motion;
        int32_t tmp = interval_max;
        interval_max = -interval_min;
        interval_min = -tmp;
    }
    // if the interval ends before the motion even begins, then we're never in
    // the interval.
    if (interval_max < 0)
        return false;
    // find the actual entry and exit periods.
    *entry_period = (interval_min + motion - 1) / motion;
    *exit_period = interval_max / motion;
    // return whether the range is empty (false) or not (true).
    return *entry_period <= *exit_period;
}

static bool position_may_eventually_be_visible_to_solution(struct solution *solution, struct vector original, struct vector current)
{
    int32_t u_entry, u_exit, v_entry, v_exit;
    if (!range_of_periods_where_motion_is_in_interval(current.u - original.u, solution->min_visible_u - original.u, solution->max_visible_u - original.u, &u_entry, &u_exit))
        return false;
    if (!range_of_periods_where_motion_is_in_interval(current.v - original.v, solution->min_visible_v - original.v, solution->max_visible_v - original.v, &v_entry, &v_exit))
        return false;
    // return true if the ranges overlap and false otherwise.
    return v_entry <= u_exit && u_entry <= v_exit;
}

static int compare_movements_by_number_of_atoms(const void *a, const void *b)
{
    return (int)((const struct movement *)a)->number_of_atoms - (int)((const struct movement *)b)->number_of_atoms;
}

static void perform_arm_instructions(struct solution *solution, struct board *board)
{
    if (solution->tape_period == 0)
        return;
    board->moving_atoms.length = 0;
    board->moving_atoms.cursor = 0;
    board->movements.length = 0;
    board->movements.cursor = 0;
    uint32_t n = solution->number_of_arms;
    for (uint32_t i = 0; i < n; ++i) {
        struct mechanism *m = &solution->arms[i];
        m->movement = zero_vector;
        m->type &= ~MOVED_GRABBED_ATOMS;
        if (board->cycle < (uint64_t)solution->arm_tape_start_cycle[i])
            continue;
        size_t index = board->cycle - (uint64_t)solution->arm_tape_start_cycle[i];
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
        struct vector track_motion = {0, 0};
        // kind of cheesy but it works!
        int32_t arm_length = abs(m->direction_u.u) | abs(m->direction_u.v);
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
            report_collision(board, m->position, "trying to extend/retract a non-piston arm");
            continue;
        } else if (inst == 'w') {
            // don't extend pistons past 3 hexes of length.
            if (arm_length == 3)
                continue;
        } else if (inst == 's') {
            // don't retract pistons below 1 hex of length.
            if (arm_length == 1)
                continue;
        } else if (inst == 't' || inst == 'g') {
            uint32_t index;
            if (!lookup_track(solution, m->position, &index)) {
                report_collision(board, m->position, "trying to move an arm along a track that isn't on a track");
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
        m->type |= MOVED_GRABBED_ATOMS;
        // next, apply the instruction to any grabbed atoms. atom movements are
        // added to a list (board->movements) and deferred until later.
        int step = angular_distance_between_grabbers(m->type);
        for (int direction = 0; direction < 6; direction += step) {
            struct vector offset = u_offset_for_direction(direction);
            struct vector p = mechanism_relative_position(*m, offset.u, offset.v, 1);
            if (inst == 'a')
                record_swing_area(board, m->position, p, 1);
            if (inst == 'd')
                record_swing_area(board, m->position, p, -1);
            struct atom_ref_at_position a = get_atom(board, *m, offset.u, offset.v);
            if (!*a.atom)
                continue;
            if (inst == 'r' && !(*a.atom & REMOVED)) {
                for (int i = 0; i < NUMBER_OF_ATOM_TYPES; ++i) {
                    if (*a.atom & ATOM_OF_TYPE(i)) {
                        board->atom_grabs[i]++;
                        break;
                    }
                }
                atom grabs = (*a.atom & GRABBED) / GRABBED_ONCE;
                *a.atom &= ~GRABBED;
                *a.atom |= (grabs + 1) * GRABBED_ONCE;
                m->type |= GRABBING_LOW_BIT << direction;
                continue;
            }
            if (!(m->type & (GRABBING_LOW_BIT << direction)) || (!(*a.atom & REMOVED) && !(*a.atom & GRABBED)))
                continue;
            if (inst == 'f' && !(*a.atom & REMOVED) && (*a.atom & GRABBED)) {
                atom grabs = (*a.atom & GRABBED) / GRABBED_ONCE;
                *a.atom &= ~GRABBED;
                *a.atom |= (grabs - 1) * GRABBED_ONCE;
                // allow conduits to transport this atom (and any atoms bonded to it).
                *a.atom |= BEING_DROPPED;
                schedule_flag_reset_if_needed(board, a.atom);
                continue;
            }
            struct movement movement = {
                .type = (m->type & PISTON) ? IS_PISTON : 0,
                .base = m->position,
                .base_rotation = m->arm_rotation,
                .grabber_offset = offset,
                .absolute_grab_position = a.position,
            };
            movement.grabber_offset.u *= arm_length;
            movement.grabber_offset.v *= arm_length;
            switch (inst) {
            case 'q': // pivot ccw
                movement.type |= PIVOT_MOVEMENT;
                movement.rotation = 1;
                move_atoms(board, a.atom, movement);
                break;
            case 'e': // pivot cw
                movement.type |= PIVOT_MOVEMENT;
                movement.rotation = -1;
                move_atoms(board, a.atom, movement);
                break;
            case 'a': // rotate ccw
                movement.type |= SWING_MOVEMENT;
                movement.rotation = 1;
                move_atoms(board, a.atom, movement);
                break;
            case 'd': // rotate cw
                movement.type |= SWING_MOVEMENT;
                movement.rotation = -1;
                move_atoms(board, a.atom, movement);
                break;
            case 'w': // extend piston
                movement.type |= PISTON_MOVEMENT;
                movement.piston_extension = 1;
                movement.translation = normalize_axis(m->direction_u);
                move_atoms(board, a.atom, movement);
                break;
            case 's': // retract piston
                movement.type |= PISTON_MOVEMENT;
                movement.piston_extension = -1;
                movement.translation = normalize_axis(m->direction_u);
                movement.translation.u = -movement.translation.u;
                movement.translation.v = -movement.translation.v;
                move_atoms(board, a.atom, movement);
                break;
            case 't': // move along track, - direction
            case 'g': // move along track, + direction
                movement.type |= TRACK_MOVEMENT;
                movement.translation = track_motion;
                move_atoms(board, a.atom, movement);
                break;
            default:
                break;
            }
        }
        // finally, transform the arm itself.
        switch (inst) {
        case 'q': // pivot ccw
        case 'e': // pivot cw
            m->pivot_parity = !m->pivot_parity;
            break;
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
            m->arm_rotation++;
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
            m->arm_rotation--;
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
            m->movement = track_motion;
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
        // the weird (uint32_t)-(int64_t) thing is to handle the INT_MIN case.
        if (m->arm_rotation > 0 && (uint32_t)m->arm_rotation > solution->maximum_absolute_arm_rotation)
            solution->maximum_absolute_arm_rotation = (uint32_t)m->arm_rotation;
        else if (m->arm_rotation < 0 && (uint32_t)-(int64_t)m->arm_rotation > solution->maximum_absolute_arm_rotation)
            solution->maximum_absolute_arm_rotation = (uint32_t)-(int64_t)m->arm_rotation;
    }
    // carry out deferred movements.
    if (board->half_cycle == 2) {
        // report a collision if any static arm is grabbing a moving atom
        for (uint32_t i = 0; i < n; ++i) {
            struct mechanism *m = &solution->arms[i];
            if (board->cycle < (uint64_t)solution->arm_tape_start_cycle[i])
                continue;
            // check whether the arm is static on this cycle.
            if (m->type & MOVED_GRABBED_ATOMS)
                continue;
            int step = angular_distance_between_grabbers(m->type);
            for (int direction = 0; direction < 6; direction += step) {
                if (!(m->type & (GRABBING_LOW_BIT << direction)))
                    continue;
                struct vector offset = u_offset_for_direction(direction);
                struct vector pos = mechanism_relative_position(*m, offset.u, offset.v, 1);
                atom *a = lookup_atom(board, pos);
                if ((*a & VALID) && (*a & REMOVED) && (*a & MOVED))
                    report_collision(board, pos, "atom moved while being held stationary by another arm");
            }
        }

        size_t atom_index = 0;

        // sort movements by number of atoms to make the collision bounding box
        // optimization work better.  this has to be done before fixing up the
        // chain atom list pointers.
        qsort(board->movements.movements, board->movements.length, sizeof(struct movement), compare_movements_by_number_of_atoms);
        size_t offset = board->moving_atoms.length;
        if (offset * 2 > board->moving_atoms.capacity) {
            board->moving_atoms.capacity = offset * 2;
            board->moving_atoms.atoms_at_positions = realloc(board->moving_atoms.atoms_at_positions, sizeof(struct atom_at_position) * board->moving_atoms.capacity);
        }
        memcpy(board->moving_atoms.atoms_at_positions + offset, board->moving_atoms.atoms_at_positions, sizeof(struct atom_at_position) * offset);
        for (size_t i = 0; i < board->movements.length; ++i) {
            struct movement *m = &board->movements.movements[i];
            memcpy(board->moving_atoms.atoms_at_positions + atom_index, board->moving_atoms.atoms_at_positions + offset + m->first_atom_index, m->number_of_atoms * sizeof(struct atom_at_position));
            m->first_atom_index = atom_index;
            atom_index += m->number_of_atoms;
        }

        int32_t maximum_rotation_distance = 1;
        atom_index = 0;
        // this is kind of terrible.  we have to fix up the movement structs to
        // look like the game is expecting (instead of the initial state, before
        // the movement, the game expects to see the final state, after the
        // movement already happened).  so, for now, we'll do two passes.  this
        // code should be cleaned up at some point.
        for (size_t i = 0; i < board->movements.length; ++i) {
            struct movement *m = &board->movements.movements[i];
            if (m->first_chain_atom != UINT32_MAX) {
                // m->prev_in_list is pointing to stack data right now; fix it so it
                // properly points to the list in the heap.
                board->chain_atoms[m->first_chain_atom].prev_in_list = &m->first_chain_atom;
            }
            // printf("movement: %d %d by %d / %d %d around %d %d (%d)\n", m.absolute_grab_position.u, m.absolute_grab_position.v, m.rotation, m.translation.u, m.translation.v, m.base.u, m.base.v, m.type);
            struct vector base = ((m->type & 3) == SWING_MOVEMENT) ? m->base : m->absolute_grab_position;
            struct vector u = u_offset_for_direction(m->rotation);
            struct vector v = v_offset_for_direction(m->rotation);
            for (size_t j = 0; j < m->number_of_atoms; ++j) {
                struct atom_at_position *ap = &board->moving_atoms.atoms_at_positions[atom_index++];
                apply_movement_to_position(base, m->translation, u, v, &ap->position);
                if (m->rotation != 0) {
                    struct vector delta = ap->position;
                    delta.u -= base.u;
                    delta.v -= base.v;
                    int32_t distance = (abs(delta.u) + abs(delta.v) + abs(delta.u + delta.v)) / 2;
                    if (distance > maximum_rotation_distance)
                        maximum_rotation_distance = distance;
                    rotate_bonds(&ap->atom, m->rotation);
                }
            }
            uint32_t chain = m->first_chain_atom;
            while (chain != UINT32_MAX) {
                struct chain_atom *ca = &board->chain_atoms[chain];
                apply_movement_to_position(base, m->translation, u, v, &ca->current_position);
                if (board->chain_mode == EXTEND_CHAIN)
                    apply_movement_to_position(base, m->translation, u, v, &ca->original_position);
                int rotation = ca->flags & CHAIN_ATOM_ROTATION;
                ca->flags &= ~CHAIN_ATOM_ROTATION;
                ca->flags |= normalize_direction(rotation + m->rotation);
                if (m->rotation > 0) {
                    uint32_t sextants = ((ca->flags & CHAIN_ATOM_SWING_SEXTANTS) >> 1) & CHAIN_ATOM_SWING_SEXTANTS;
                    sextants |= 1u << (CHAIN_ATOM_SWING_SEXTANTS_SHIFT + 5);
                    ca->flags &= ~CHAIN_ATOM_SWING_SEXTANTS;
                    ca->flags |= sextants;
                } else if (m->rotation < 0) {
                    uint32_t sextants = ((ca->flags & CHAIN_ATOM_SWING_SEXTANTS) << 1) & CHAIN_ATOM_SWING_SEXTANTS;
                    sextants |= 1u << CHAIN_ATOM_SWING_SEXTANTS_SHIFT;
                    ca->flags &= ~CHAIN_ATOM_SWING_SEXTANTS;
                    ca->flags |= sextants;
                }
                chain = ca->next_in_list;
            }
            if ((m->type & 3) == TRACK_MOVEMENT) {
                m->base.u += m->translation.u;
                m->base.v += m->translation.v;
            } else if ((m->type & 3) == PISTON_MOVEMENT) {
                struct vector g = normalize_axis(m->grabber_offset);
                m->grabber_offset.u = g.u ? g.u * (m->grabber_offset.u / g.u + m->piston_extension) : 0;
                m->grabber_offset.v = g.v ? g.v * (m->grabber_offset.v / g.v + m->piston_extension) : 0;
            }
            apply_movement_to_position(base, m->translation, u, v, &m->absolute_grab_position);
            if ((m->type & 3) == SWING_MOVEMENT)
                m->base_rotation += m->rotation;
        }
        double collision_increment = 0.25 / pow(2, round(log(maximum_rotation_distance) / log(2)));
        if (!(collision_increment <= 0.125))
            collision_increment = 0.125;
        struct vector collision_location;
        if (collision(solution, board, (float)collision_increment, &collision_location))
            report_collision(board, collision_location, "collision during motion phase");
        if (board->collision_check_limit > 0 && board->collision_checks > board->collision_check_limit)
            report_collision(board, zero_vector, "solution reached limit for maximum number of collision checks");
        atom_index = 0;
        for (size_t i = 0; i < board->movements.length; ++i) {
            struct movement m = board->movements.movements[i];
            for (size_t j = 0; j < m.number_of_atoms; ++j) {
                struct atom_at_position ap = board->moving_atoms.atoms_at_positions[atom_index++];
                if (position_may_be_visible_to_solution(solution, ap.position))
                    ap.atom &= ~IS_CHAIN_ATOM;
                insert_atom(board, ap.position, VALID | ap.atom, "atom moved on top of another atom");
            }
            uint32_t chain = m.first_chain_atom;
            while (chain != UINT32_MAX) {
                struct chain_atom ca = board->chain_atoms[chain];
                if (position_may_be_visible_to_solution(solution, ca.current_position)) {
                    if (board->chain_mode == EXTEND_CHAIN)
                        board->chain_will_become_visible = true;
                    move_chain_atom_to_list(board, chain, 0);
                } else if (board->chain_mode == EXTEND_CHAIN && position_may_eventually_be_visible_to_solution(solution, ca.original_position, ca.current_position)) {
                    board->chain_will_become_visible = true;
                    *lookup_atom_in_grid(&board->grid, ca.current_position) &= ~IS_CHAIN_ATOM;
                    move_chain_atom_to_list(board, chain, 0);
                } else
                    add_chain_atom_to_table(board, chain);
                chain = ca.next_in_list;
            }
        }
    }
}

static void mark_arm_area(struct solution *solution, struct board *board)
{
    // xx definitely do this in a cleaner way...
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        struct mechanism *m = &solution->arms[i];
        int step = angular_distance_between_grabbers(m->type);
        for (int direction = 0; direction < 6; direction += step) {
            struct vector offset = u_offset_for_direction(direction);
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
    int rotation = direction_for_offset(m.direction_u);
    while (*atom_cursor < *atom_next) {
        struct atom_at_position ap = conduit->atoms[*atom_cursor];
        for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
            if (!(ap.atom & (BOND_LOW_BITS << normalize_direction(bond_direction + rotation)) & ~RECENT_BONDS))
                continue;
            struct vector p = ap.position;
            struct vector d = u_offset_for_direction(bond_direction);
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
            if ((*a & VALID) && !(*a & REMOVED) && !(*a & BEING_PRODUCED) && (*a & BEING_DROPPED) && !(*a & GRABBED)) {
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
                    // reset the visited flag and the cursor.
                    for (uint32_t k = saved_cursor; k < atom_next; ++k) {
                        struct vector p = conduit->atoms[k].position;
                        p = mechanism_relative_position(m, p.u, p.v, 1);
                        *lookup_atom(board, p) &= ~VISITED;
                    }
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
        for (uint32_t j = 0; j < n; ++j) {
            atom input = io->atoms[j].atom;
            struct vector p = io->atoms[j].position;
            insert_atom(board, p, VALID | input, "overlapped inputs");
        }
    }
    board->active_input_or_output = UINT32_MAX;
}

static int32_t gcd(int32_t a, int32_t b)
{
    while (b != 0) {
        int32_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static void match_repeating_output_with_chain_atoms(struct board *board, struct input_output *io)
{
    struct atom_at_position placeholder = io->original_atoms[io->number_of_original_atoms - 1];
    int32_t offset = placeholder.position.u - io->repetition_origin.u;
    if (offset <= 0 || placeholder.position.v != io->repetition_origin.v)
        return;
    int32_t maximum_feed_rate = 0;
    for (uint32_t i = 0; i < io->number_of_original_atoms - 1; ++i) {
        struct atom_at_position output = io->original_atoms[i];
        if (vectors_equal(output.position, io->repetition_origin))
            output.atom |= placeholder.atom & ALL_BONDS;
        uint32_t repeats_after = UINT32_MAX;
        for (uint32_t j = REPEATING_OUTPUT_REPETITIONS; j < repeats_after; ++j) {
            struct vector p = output.position;
            p.u += j * offset;
            atom a = *lookup_atom(board, p);
            if (!(a & VALID) || (a & REMOVED) || (a & BEING_PRODUCED) || (a & IS_CHAIN_ATOM)) {
                // look for a repeating chain atom that will fill in this hex.
                bool found_chain_atom = false;
                for (uint32_t k = 0; k < board->number_of_chain_atoms; ++k) {
                    struct chain_atom ca = board->chain_atoms[k];
                    if (!ca.prev_in_list || !(ca.flags & CHAIN_ATOM_IN_REPEATING_SEGMENT))
                        continue;
                    if (ca.current_position.v != p.v || ca.current_position.u > p.u)
                        continue;
                    int32_t chain_offset = ca.current_position.u - ca.original_position.u;
                    if (chain_offset <= 0 || ca.current_position.v != ca.original_position.v)
                        continue;
                    if ((p.u - ca.current_position.u) % chain_offset != 0)
                        continue;
                    // we only have to check up to the offset where this chain
                    // atom appears in the same spot in the output.
                    if (repeats_after == UINT32_MAX)
                        repeats_after = j + chain_offset / gcd(offset, chain_offset);
                    if (chain_offset > maximum_feed_rate)
                        maximum_feed_rate = chain_offset;
                    found_chain_atom = true;
                    a = *lookup_atom(board, ca.current_position);
                    break;
                }
                if (!found_chain_atom)
                    return;
            }
            if ((a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS))) != output.atom)
                return;
        }
    }
    // the entire product matched; update the maximum feed rate.
    if (maximum_feed_rate > io->maximum_feed_rate)
        io->maximum_feed_rate = maximum_feed_rate;
}

static bool mark_output_position(struct board *board, struct vector pos)
{
    atom *a = lookup_atom(board, pos);
    if ((*a & VISITED) || !(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED))
        return false;
    *a |= VISITED;
    size_t capacity = board->marked.capacity;
    while (board->marked.length >= capacity)
        capacity = 4 * (capacity + 12) / 3;
    if (capacity != board->marked.capacity) {
        struct vector *positions = realloc(board->marked.positions,
         sizeof(struct movement) * capacity);
        if (!positions)
            abort();
        board->marked.positions = positions;
        board->marked.capacity = capacity;
    }
    board->marked.positions[board->marked.length++] = pos;
    return true;
}

static void consume_outputs(struct solution *solution, struct board *board)
{
    size_t active = board->active_input_or_output;
    size_t i = active == UINT32_MAX ? 0 : active;
    for (; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        if (!(io->type & OUTPUT))
            continue;
        bool fails_on_wrong_output = board->fails_on_wrong_output_mask & (1ULL << io->puzzle_index);
        bool wrong_output = fails_on_wrong_output;
        board->marked.length = 0;
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
                wrong_output = false;
                break;
            }
            // variable outputs match any atom.
            if (output & VARIABLE_OUTPUT) {
                output &= ~VARIABLE_OUTPUT;
                output |= *a & ANY_ATOM;
            }
            if (io->type & REPEATING_OUTPUT) {
                uint64_t bond_mask = (((output & RECENT_BONDS) >> RECENT_BOND) * BOND_LOW_BITS) & ~RECENT_BONDS;
                if ((*a & (ANY_ATOM | bond_mask)) != (output & (ANY_ATOM | bond_mask))) {
                    // printf("did not match at %d %d: %llx vs %llx\n", io->atoms[j].position.u, io->atoms[j].position.v, (*a & (ANY_ATOM | bond_mask)), output & (ANY_ATOM | bond_mask));
                    match = false;
                    break;
                }
                mark_output_position(board, io->atoms[j].position);
            } else if (fails_on_wrong_output) {
                if ((*a & (ALL_BONDS & ~RECENT_BONDS)) != (output & ALL_BONDS & ~RECENT_BONDS) || (*a & GRABBED)) {
                    match = false;
                    wrong_output = false;
                    break;
                }
                if ((*a & ANY_ATOM) != (output & ANY_ATOM)) {
                    match = false;
                    // this could be a wrong output.
                }
            } else {
                if ((*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS))) != output || (*a & GRABBED)) {
                    // printf("did not match at %d %d: %llx vs %llx\n", io->atoms[j].position.u, io->atoms[j].position.v, (*a & (ANY_ATOM | (ALL_BONDS & ~RECENT_BONDS))), output);
                    match = false;
                    break;
                }
            }
            // printf("match!\n");
        }
        bool wrong_output_bonds = false;
        if (board->fails_on_wrong_output_bonds_mask & (1ULL << io->puzzle_index)) {
            wrong_output_bonds = true;
            for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
                if (io->atoms[j].atom & REPEATING_OUTPUT_PLACEHOLDER)
                    continue;
                atom *a = lookup_atom(board, io->atoms[j].position);
                if (!(*a & VALID) || (*a & REMOVED) || (*a & BEING_PRODUCED) || (*a & GRABBED)) {
                    wrong_output_bonds = false;
                    break;
                }
                mark_output_position(board, io->atoms[j].position);
            }
            // if wrong_output_bonds is set at this point, then the molecule
            // covers the output completely.  the goal now is to determine
            // whether there are any bonds leading out of the output shape.  if
            // so, the footprint doesn't match, so the output doesn't have wrong
            // bonds.
            for (uint32_t j = 0; !match && wrong_output_bonds && j < io->number_of_atoms; ++j) {
                if (io->atoms[j].atom & REPEATING_OUTPUT_PLACEHOLDER)
                    continue;
                atom a = *lookup_atom(board, io->atoms[j].position);
                for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
                    if (!(a & (BOND_LOW_BITS << bond_direction) & ~RECENT_BONDS))
                        continue;
                    struct vector next = io->atoms[j].position;
                    struct vector d = u_offset_for_direction(bond_direction);
                    next.u += d.u;
                    next.v += d.v;
                    atom b = *lookup_atom(board, next);
                    if (!(b & VISITED)) {
                        wrong_output_bonds = false;
                        break;
                    }
                }
            }
        }
        // validating infinite products requires visiting all the atoms in the
        // molecule.  that's what this loop does.
        if (io->type & REPEATING_OUTPUT) {
            size_t cursor = 0;
            while (match && cursor < board->marked.length) {
                struct vector p = board->marked.positions[cursor++];
                atom a = *lookup_atom(board, p);
                if (a & IS_CHAIN_ATOM && board->chain_mode == EXTEND_CHAIN) {
                    // ensure atoms stay within the vertical region forever.
                    struct chain_atom ca = board->chain_atoms[lookup_chain_atom(board, p)];
                    if (ca.current_position.v - ca.original_position.v != 0) {
                        match = false;
                        break;
                    }
                }
                for (int bond_direction = 0; bond_direction < 6; ++bond_direction) {
                    if (!(a & (BOND_LOW_BITS << bond_direction) & ~RECENT_BONDS))
                        continue;
                    struct vector next = p;
                    struct vector d = u_offset_for_direction(bond_direction);
                    next.u += d.u;
                    next.v += d.v;
                    if (!mark_output_position(board, next))
                        continue;
                    // atoms cannot appear outside the vertical bounds of the
                    // infinite product.
                    if (next.v < io->min_v || next.v > io->max_v) {
                        match = false;
                        break;
                    }
                    // atoms cannot appear between atoms in a row of the infinite
                    // product.  the row_min_v and row_max_v arrays track the
                    // disallowed range of positions for each row.
                    size_t row = next.v - io->min_v;
                    if (next.u >= io->row_min_u[row] && next.u <= io->row_max_u[row]) {
                        match = false;
                        break;
                    }
                }
            }
        }
        // reset the visited flag for all the marked atoms.
        for (size_t j = 0; j < board->marked.length; ++j) {
            atom *a = lookup_atom(board, board->marked.positions[j]);
            *a &= ~VISITED;
        }
        if (match && (io->type & REPEATING_OUTPUT) && board->chain_mode == EXTEND_CHAIN)
            match_repeating_output_with_chain_atoms(board, io);
        // if the output is a match, first trigger an interrupt if necessary.
        // then remove the output and increment the output counter.
        if (match) {
            if (i != active && (io->type & INTERRUPT)) {
                board->active_input_or_output = i;
                return;
            }
            if (io->type & REPEATING_OUTPUT)
                io->number_of_outputs = io->number_of_repetitions * io->outputs_per_repetition;
            else {
                for (uint32_t j = 0; j < io->number_of_atoms; ++j)
                    remove_atom(board, (struct atom_ref_at_position){ lookup_atom(board, io->atoms[j].position), io->atoms[j].position });
                io->number_of_outputs++;
            }
        } else if (wrong_output) {
            report_collision(board, io->atoms[0].position, "output didn't match");
            board->wrong_output_index = i;
        } else if (wrong_output_bonds) {
            report_collision(board, io->atoms[0].position, "output bonds didn't match");
            board->wrong_output_index = i;
        }
    }
    board->active_input_or_output = UINT32_MAX;
}

static void reset_temporary_flags(struct board *board)
{
    for (uint32_t i = 0; i < board->flag_reset_length; ++i)
        *board->flag_reset[i] &= ~TEMPORARY_FLAGS;
    board->flag_reset_length = 0;
    for (uint32_t i = 0; i < board->number_of_overlapped_atoms; ++i)
        board->overlapped_atoms[i].atom &= ~TEMPORARY_FLAGS;
}

static void check_completion(struct solution *solution, struct board *board)
{
    uint64_t min = UINT64_MAX;
    for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        if (!(solution->inputs_and_outputs[i].type & OUTPUT))
            continue;
        uint64_t count = solution->inputs_and_outputs[i].number_of_outputs;
        if (count < min)
            min = count;
    }
    board->complete = min >= solution->target_number_of_outputs;
    if (min > board->output_cycles_capacity) {
        board->output_cycles_capacity = (min * 3) / 2;
        if (board->output_cycles_capacity < 6)
            board->output_cycles_capacity = 6;
        board->output_cycles = realloc(board->output_cycles, board->output_cycles_capacity * sizeof(uint64_t));
    }
    for (uint64_t i = board->number_of_output_cycles; i < min; ++i)
        board->output_cycles[i] = board->cycle;
    board->number_of_output_cycles = min;
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
        if (board->half_cycle == 2 && board->number_of_overlapped_atoms > 0)
            report_collision(board, board->overlapped_atoms[0].position, "overlapping atoms before motion phase");
        perform_arm_instructions(solution, board);
        if (board->half_cycle == 2 && board->number_of_overlapped_atoms > 0)
            report_collision(board, board->overlapped_atoms[0].position, "overlapping atoms after motion phase");
        mark_arm_area(solution, board);
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
    board->cycle++;
    board->half_cycle = 1;
    check_completion(solution, board);
    return FINISHED_CYCLE;
}

bool repeat_molecule(struct input_output *io, uint32_t repetitions, const char **error)
{
    if (io->number_of_original_atoms == 0) {
        *error = "puzzle contains an empty infinite product";
        return false;
    }
    struct atom_at_position placeholder = io->original_atoms[io->number_of_original_atoms - 1];
    struct vector offset = placeholder.position;
    offset.u -= io->repetition_origin.u;
    offset.v -= io->repetition_origin.v;
    placeholder.position.u += offset.u * (repetitions - 1);
    placeholder.position.v += offset.v * (repetitions - 1);
    struct atom_at_position *atoms = calloc((io->number_of_original_atoms - 1) * repetitions + 1,
     sizeof(io->atoms[0]));
    for (uint32_t i = 0; i < repetitions; ++i) {
        for (uint32_t j = 0; j < io->number_of_original_atoms - 1; ++j) {
            struct atom_at_position *a = &atoms[(io->number_of_original_atoms - 1) * i + j];
            *a = io->original_atoms[j];
            bool origin = vectors_equal(a->position, io->repetition_origin);
            a->position.u += i * offset.u;
            a->position.v += i * offset.v;
            a->atom = io->original_atoms[j].atom;
            // add the incoming bonds from the previous monomer.
            if (origin && i > 0)
                a->atom |= placeholder.atom & ALL_BONDS;
        }
    }
    atoms[(io->number_of_original_atoms - 1) * repetitions] = placeholder;
    free(io->atoms);
    io->atoms = atoms;
    io->number_of_atoms = (io->number_of_original_atoms - 1) * repetitions + 1;
    io->number_of_repetitions = repetitions;

    io->min_v = INT32_MAX;
    io->max_v = INT32_MIN;
    for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
        if (io->atoms[j].atom & REPEATING_OUTPUT_PLACEHOLDER)
            continue;
        struct vector p = io->atoms[j].position;
        if (p.v < io->min_v)
            io->min_v = p.v;
        if (p.v > io->max_v)
            io->max_v = p.v;
    }
    if ((int64_t)io->max_v - (int64_t)io->min_v > 99999) {
        *error = "solution contains an infinite product with too many rows";
        return false;
    }
    size_t rows = io->max_v - io->min_v + 1;
    io->row_min_u = realloc(io->row_min_u, rows * sizeof(int32_t));
    io->row_max_u = realloc(io->row_max_u, rows * sizeof(int32_t));
    for (size_t j = 0; j < rows; ++j) {
        io->row_min_u[j] = INT32_MAX;
        io->row_max_u[j] = INT32_MIN;
    }
    for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
        if (io->atoms[j].atom & REPEATING_OUTPUT_PLACEHOLDER)
            continue;
        struct vector p = io->atoms[j].position;
        size_t row = p.v - io->min_v;
        if (p.u < io->row_min_u[row])
            io->row_min_u[row] = p.u;
        if (p.u > io->row_max_u[row])
            io->row_max_u[row] = p.u;
    }
    for (uint32_t j = 0; j < io->number_of_atoms; ++j) {
        if (io->atoms[j].atom & REPEATING_OUTPUT_PLACEHOLDER)
            continue;
        struct vector p = io->atoms[j].position;
        size_t row = p.v - io->min_v;
        if (p.u + 1 <= io->row_max_u[row])
            io->atoms[j].atom |= 1ULL << (RECENT_BOND + 0);
        if (p.v != io->max_v) {
            size_t neighbor = p.v + 1 - io->min_v;
            if (p.u >= io->row_min_u[neighbor] && p.u <= io->row_max_u[neighbor])
                io->atoms[j].atom |= 1ULL << (RECENT_BOND + 1);
            if (p.u - 1 >= io->row_min_u[neighbor] && p.u - 1 <= io->row_max_u[neighbor])
                io->atoms[j].atom |= 1ULL << (RECENT_BOND + 2);
        }
        if (p.u - 1 >= io->row_min_u[row])
            io->atoms[j].atom |= 1ULL << (RECENT_BOND + 3);
        if (p.v != io->min_v) {
            size_t neighbor = p.v - 1 - io->min_v;
            if (p.u >= io->row_min_u[neighbor] && p.u <= io->row_max_u[neighbor])
                io->atoms[j].atom |= 1ULL << (RECENT_BOND + 4);
            if (p.u + 1 >= io->row_min_u[neighbor] && p.u + 1 <= io->row_max_u[neighbor])
                io->atoms[j].atom |= 1ULL << (RECENT_BOND + 5);
        }
    }
    return true;
}

struct footprint {
    enum mechanism_type type;
    const struct vector *hexes;
};

static struct footprint footprints[] = {
    {
        .type = CALCIFICATION,
        .hexes = (const struct vector[]){
            { 0, 0 },
        }
    },
    {
        .type = ANIMISMUS,
        .hexes = (const struct vector[]){
            { 0, 1 },
            { 1, 0 },
            { 1, -1 },
            { 0, 0 },
        }
    },
    {
        .type = PROJECTION,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 0 },
        }
    },
    {
        .type = DISPERSION,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 1, -1 },
            { 0, -1 },
            { -1, 0 },
            { 0, 0 },
        }
    },
    {
        .type = PURIFICATION,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 1 },
            { 0, 0 },
        }
    },
    {
        .type = DUPLICATION,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 0 },
        }
    },
    {
        .type = UNIFICATION,
        .hexes = (const struct vector[]){
            { 0, 1 },
            { -1, 1 },
            { 0, -1 },
            { 1, -1 },
            { 0, 0 },
        }
    },
    {
        .type = BONDING,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 0 },
        }
    },
    {
        .type = UNBONDING,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 0 },
        }
    },
    {
        .type = TRIPLEX_BONDING,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 1 },
            { 0, 0 },
        }
    },
    {
        .type = MULTI_BONDING,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, -1 },
            { -1, 1 },
            { 0, 0 },
        }
    },
    {
        .type = DISPOSAL,
        .hexes = (const struct vector[]){
            { 1, 0 },
            { 0, 1 },
            { -1, 1 },
            { -1, 0 },
            { 0, -1 },
            { 1, -1 },
            { 0, 0 },
        }
    },
    {
        .type = EQUILIBRIUM,
        .hexes = (const struct vector[]){
            { 0, 0 },
        }
    },
};

const struct vector *glyph_footprint(enum mechanism_type type)
{
    for (int i = 0; i < sizeof(footprints)/sizeof(footprints[0]); ++i) {
        if (type & footprints[i].type)
            return footprints[i].hexes;
    }
    static const struct vector trivial_footprint[] = { { 0, 0 } };
    return trivial_footprint;
}

static void create_van_berlo_atom(struct board *board, struct mechanism m, int32_t du, int32_t dv, atom element)
{
    struct vector p = mechanism_relative_position(m, du, dv, 1);
    insert_atom(board, p, VALID | GRABBED_ONCE | VAN_BERLO_ATOM | element, "van berlo overlap");
}

void initial_setup(struct solution *solution, struct board *board, uint32_t initial_board_size)
{
    board->overlap = solution->track_self_overlap;
    board->half_cycle = 1;
    board->active_input_or_output = UINT32_MAX;
    board->wrong_output_index = SIZE_MAX;
    // initialize the board array / hash table.
    // this division is just to keep the hash table size from starting too big
    // and ruining memory locality.
    uint32_t initial_hashtable_size = initial_board_size / 8;
    if (initial_hashtable_size < 1)
        initial_hashtable_size = 1;
    if (initial_hashtable_size > 99999)
        initial_hashtable_size = 99999;
    board->grid.board = board;
    rehash(&board->grid, initial_hashtable_size);
    // place cabinet walls first; wall-wall overlap isn't counted.
    for (uint32_t i = 0; i < solution->number_of_cabinet_walls; ++i)
        mark_used_area(board, solution->cabinet_walls[i]);
    for (uint32_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        for (uint32_t j = 0; j < io->number_of_atoms; ++j)
            mark_used_area_with_overlap(board, io->atoms[j].position, &board->overlap);
    }
    for (uint32_t i = 0; i < solution->number_of_glyphs; ++i) {
        struct mechanism m = solution->glyphs[i];
        if (m.type & CONDUIT)
            continue;
        const struct vector *footprint = glyph_footprint(m.type);
        for (int j = 0; ; ++j) {
            struct vector p = footprint[j];
            mark_used_area_with_overlap(board, mechanism_relative_position(m, p.u, p.v, 1), &board->overlap);
            if (vectors_equal(p, zero_vector))
                break;
        }
    }
    for (uint32_t i = 0; i < solution->number_of_conduits; ++i) {
        struct conduit *conduit = &solution->conduits[i];
        struct mechanism m = solution->glyphs[conduit->glyph_index];
        for (uint32_t j = 0; j < conduit->number_of_positions; ++j) {
            struct vector p = conduit->positions[j];
            mark_used_area_with_overlap(board, mechanism_relative_position(m, p.u, p.v, 1), &board->overlap);
        }
    }
    for (size_t i = 0; i < solution->number_of_arms; ++i)
        mark_used_area_with_overlap(board, solution->arms[i].position, &board->overlap);
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        // use the VISITED flag to mark arm base positions that track is allowed to occupy.
        *lookup_atom(board, solution->arms[i].position) |= VISITED;
    }
    for (uint32_t i = 0; i < solution->track_table_size; ++i) {
        struct vector p = solution->track_positions[i];
        if (p.u == INT32_MIN && p.v == INT32_MIN)
            continue;
        mark_used_area_with_overlap(board, p, &board->overlap);
    }
    for (size_t i = 0; i < solution->number_of_arms; ++i)
        *lookup_atom(board, solution->arms[i].position) &= ~VISITED;
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        if (!(solution->arms[i].type & VAN_BERLO))
            continue;
        solution->arms[i].type |= GRABBING;
        solution->arms[i].type |= GRABBING_EVERYTHING;
        create_van_berlo_atom(board, solution->arms[i], 1, 0, SALT);
        create_van_berlo_atom(board, solution->arms[i], 0, 1, WATER);
        create_van_berlo_atom(board, solution->arms[i], -1, 1, AIR);
        create_van_berlo_atom(board, solution->arms[i], -1, 0, SALT);
        create_van_berlo_atom(board, solution->arms[i], 0, -1, FIRE);
        create_van_berlo_atom(board, solution->arms[i], 1, -1, EARTH);
    }
    mark_arm_area(solution, board);
    // van berlo's wheel can block inputs.
    flag_blocked_inputs(solution, board);
    spawn_inputs(solution, board);
}

void destroy(struct solution *solution, struct board *board)
{
    if (solution) {
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
        free(solution->cabinet_walls);
        for (size_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
            if (solution->inputs_and_outputs[i].original_atoms != solution->inputs_and_outputs[i].atoms)
                free(solution->inputs_and_outputs[i].original_atoms);
            free(solution->inputs_and_outputs[i].atoms);
            free(solution->inputs_and_outputs[i].row_min_u);
            free(solution->inputs_and_outputs[i].row_max_u);
        }
        free(solution->inputs_and_outputs);
        memset(solution, 0, sizeof(*solution));
    }
    if (board) {
        free(board->grid.atoms_at_positions);
        free(board->flag_reset);
        free(board->movements.movements);
        free(board->moving_atoms.atoms_at_positions);
        free(board->marked.positions);
        free(board->overlapped_atoms);
        free(board->chain_atoms);
        free(board->chain_atom_table);
        free(board->output_cycles);
        for (uint32_t i = 0; i < board->number_of_area_directions; ++i)
            free(board->area_directions[i].footprint_at_infinity.atoms_at_positions);
        free(board->area_directions);
        memset(board, 0, sizeof(*board));
    }
}

uint32_t used_area(struct board *board)
{
    return board->area;
}

// appendix A -- hash table functions.

bool lookup_track(struct solution *solution, struct vector query, uint32_t *index)
{
    if (solution->track_table_size == 0)
        return false;
    uint32_t hash = vector_hash(query);
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

void move_chain_atom_to_list(struct board *board, uint32_t chain_atom_index, uint32_t *list)
{
    if (chain_atom_index == UINT32_MAX)
        return;
    struct chain_atom *m = &board->chain_atoms[chain_atom_index];
    if (m->next_in_list != UINT32_MAX)
        board->chain_atoms[m->next_in_list].prev_in_list = m->prev_in_list;
    if (m->prev_in_list)
        *m->prev_in_list = m->next_in_list;
    if (list) {
        if (*list == chain_atom_index)
            abort();
        m->prev_in_list = list;
        m->next_in_list = *list;
        if (*list != UINT32_MAX)
            board->chain_atoms[*list].prev_in_list = &m->next_in_list;
        *list = chain_atom_index;
    } else {
        m->prev_in_list = 0;
        m->next_in_list = UINT32_MAX;
    }
}

void add_chain_atom_to_table(struct board *board, uint32_t chain_atom_index)
{
    if (board->chain_atom_table_size == 0)
        return;
    struct vector p = board->chain_atoms[chain_atom_index].current_position;
    uint32_t hash = vector_hash(p);
    uint32_t mask = board->chain_atom_table_size - 1;
    move_chain_atom_to_list(board, chain_atom_index, &board->chain_atom_table[hash & mask]);
}

uint32_t lookup_chain_atom(struct board *board, struct vector query)
{
    if (board->chain_atom_table_size == 0)
        return UINT32_MAX;
    uint32_t hash = vector_hash(query);
    uint32_t mask = board->chain_atom_table_size - 1;
    uint32_t index = board->chain_atom_table[hash & mask];
    while (index != UINT32_MAX) {
        if (vectors_equal(query, board->chain_atoms[index].current_position))
            return index;
        index = board->chain_atoms[index].next_in_list;
    }
    return index;
}

static struct atom_at_position *lookup_atom_at_position_in_array(struct atom_grid *grid, struct vector query)
{
    return &grid->atoms_at_positions[(query.u - GRID_ARRAY_MIN) * (GRID_ARRAY_MAX - GRID_ARRAY_MIN) + query.v - GRID_ARRAY_MIN];
}

__attribute__((noinline))
static struct atom_at_position *lookup_atom_at_position_in_hashtable(struct atom_grid *grid, struct vector query)
{
    uint32_t hash = vector_hash(query);
    while (true) {
        uint32_t mask = grid->hash_capacity - 1;
        uint32_t index = hash & mask;
        uint32_t misses = 0;
        while (true) {
            struct atom_at_position *a = &grid->atoms_at_positions[GRID_ARRAY_PREFIX + index];
            if (!a->atom)
                return a;
            if (vectors_equal(a->position, query))
                return a;
            if (++misses == 6) {
                rehash(grid, grid->hash_capacity * 2);
                break;
            }
            index = (index + 1) & mask;
        }
        // if execution reaches this point, the board was just rehashed, so we
        // have to retry the lookup.
    }
}

struct atom_at_position *lookup_atom_at_position(struct atom_grid *grid, struct vector query)
{
    if (grid->hash_capacity == 0)
        rehash(grid, 1);
    if (query.u >= GRID_ARRAY_MIN && query.u < GRID_ARRAY_MAX && query.v >= GRID_ARRAY_MIN && query.v < GRID_ARRAY_MAX)
        return lookup_atom_at_position_in_array(grid, query);
    else
        return lookup_atom_at_position_in_hashtable(grid, query);
}

atom *lookup_atom(struct board *board, struct vector query)
{
    return &lookup_atom_at_position(&board->grid, query)->atom;
}

atom *lookup_atom_in_grid(struct atom_grid *grid, struct vector query)
{
    return &lookup_atom_at_position(grid, query)->atom;
}

static atom mark_used_area_with_overlap(struct board *board, struct vector point, uint64_t *overlap)
{
    struct atom_at_position *a = lookup_atom_at_position(&board->grid, point);
    if (a->atom & VALID) {
        if (overlap && !(a->atom & VISITED))
            (*overlap)++;
        return a->atom;
    }
    a->position = point;
    a->atom = VALID | REMOVED;
    board->area++;
    return a->atom;
}

atom mark_used_area(struct board *board, struct vector point)
{
    return mark_used_area_with_overlap(board, point, 0);
}

__attribute__((noinline))
static void insert_overlapped_atom(struct board *board, struct atom_at_position *existing_atom, atom atom, const char *collision_reason)
{
    size_t capacity = board->overlapped_atoms_capacity;
    while (board->number_of_overlapped_atoms >= capacity)
        capacity = 4 * (capacity + 12) / 3;
    if (capacity != board->number_of_overlapped_atoms) {
        struct atom_at_position *atoms_at_positions = realloc(board->overlapped_atoms,
         sizeof(struct atom_at_position) * capacity);
        if (!atoms_at_positions)
            abort();
        board->overlapped_atoms = atoms_at_positions;
        board->overlapped_atoms_capacity = capacity;
    }
    board->overlapped_atoms[board->number_of_overlapped_atoms++] = (struct atom_at_position){ .atom = atom, .position = existing_atom->position };
    existing_atom->atom |= OVERLAPS_ATOMS;
}

void insert_atom(struct board *board, struct vector query, atom atom, const char *collision_reason)
{
    struct atom_at_position *a = lookup_atom_at_position(&board->grid, query);
    if ((a->atom & VALID) && !(a->atom & REMOVED))
        insert_overlapped_atom(board, a, atom, collision_reason);
    else {
        a->position = query;
        if (!(a->atom & REMOVED))
            board->area++;
        a->atom = atom;
    }
}

__attribute__((noinline))
static void rehash(struct atom_grid *grid, uint32_t size)
{
    uint32_t n = grid->hash_capacity;
    if (n < 16)
        n = 16;
    while (n < size)
        n *= 2;
    if (n == grid->hash_capacity)
        return;
    struct atom_grid old = *grid;
    grid->atoms_at_positions = calloc(GRID_ARRAY_PREFIX + n, sizeof(*grid->atoms_at_positions));
    grid->hash_capacity = n;
    if (grid->board)
        grid->board->flag_reset_length = 0;
    if (!old.atoms_at_positions)
        return;
    for (uint32_t i = 0; i < GRID_ARRAY_PREFIX; ++i) {
        struct atom_at_position a = old.atoms_at_positions[i];
        if (!(a.atom & VALID))
            continue;
        struct atom_at_position *b = lookup_atom_at_position_in_array(grid, a.position);
        *b = a;
        if (grid->board)
            schedule_flag_reset_if_needed(grid->board, &b->atom);
    }
    for (uint32_t i = 0; i < old.hash_capacity; ++i) {
        struct atom_at_position a = old.atoms_at_positions[GRID_ARRAY_PREFIX + i];
        if (!(a.atom & VALID))
            continue;
        struct atom_at_position *b = lookup_atom_at_position_in_hashtable(grid, a.position);
        *b = a;
        if (grid->board)
            schedule_flag_reset_if_needed(grid->board, &b->atom);
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
