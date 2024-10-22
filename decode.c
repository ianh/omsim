#include "decode.h"

#include "parse.h"
#include "sim.h"
#include <stdio.h>
#include <stdlib.h>

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
    case 'C': abort(); // repeat and reset instructions are handled separately.
    case 'X': abort();
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

static uint64_t parts_available_bits_for_part_name(struct byte_string part_name)
{
    if (byte_string_is(part_name, "glyph-calcification"))
        return 1ull << 12;
    else if (byte_string_is(part_name, "glyph-life-and-death"))
        return 1ull << 16;
    else if (byte_string_is(part_name, "glyph-projection"))
        return 1ull << 14;
    else if (byte_string_is(part_name, "glyph-dispersion"))
        return 1ull << 18;
    else if (byte_string_is(part_name, "glyph-purification"))
        return 1ull << 15;
    else if (byte_string_is(part_name, "glyph-duplication"))
        return 1ull << 13;
    else if (byte_string_is(part_name, "glyph-unification"))
        return 1ull << 18;
    else if (byte_string_is(part_name, "bonder"))
        return 1ull << 8;
    else if (byte_string_is(part_name, "unbonder"))
        return 1ull << 9;
    else if (byte_string_is(part_name, "bonder-prisma"))
        return 1ull << 11;
    else if (byte_string_is(part_name, "bonder-speed"))
        return 1ull << 10;
    else if (byte_string_is(part_name, "glyph-disposal"))
        return 1ull << 17;
    else if (byte_string_is(part_name, "glyph-marker"))
        return 1ull << 1;
    else if (byte_string_is(part_name, "arm1"))
        return 1ull << 0;
    else if (byte_string_is(part_name, "arm2"))
        return 1ull << 1;
    else if (byte_string_is(part_name, "arm3"))
        return 1ull << 1;
    else if (byte_string_is(part_name, "arm6"))
        return 1ull << 1;
    else if (byte_string_is(part_name, "piston"))
        return 1ull << 2;
    else if (byte_string_is(part_name, "baron"))
        return 1ull << 28;
    else if (byte_string_is(part_name, "track"))
        return 1ull << 3;
    else
        return 0;
}

static uint64_t parts_available_bits_for_instruction(char inst)
{
    switch (inst) {
    case 'R': return 1ull << 22;
    case 'r': return 1ull << 22;
    case 'E': return 1ull << 2;
    case 'e': return 1ull << 2;
    case 'G': return 1ull << 23;
    case 'g': return 1ull << 22;
    case 'P': return 1ull << 26;
    case 'p': return 1ull << 26;
    case 'A': return 1ull << 3;
    case 'a': return 1ull << 3;
    case 'C': return 1ull << 25;
    case 'X': return 1ull << 24;
    case 'O': return 1ull << 25;
    default: return 0;
    }
}

static int compare_instructions_by_index(const void *aa, const void *bb)
{
    const struct solution_instruction *a = aa;
    const struct solution_instruction *b = bb;
    if (a->index < b->index)
        return -1;
    if (a->index > b->index)
        return 1;
    return 0;
}

static int compare_conduits_by_id(const void *aa, const void *bb)
{
    const struct conduit *a = aa;
    const struct conduit *b = bb;
    if (a->id < b->id)
        return -1;
    if (a->id > b->id)
        return 1;
    if (a->glyph_index < b->glyph_index)
        return -1;
    if (a->glyph_index > b->glyph_index)
        return 1;
    return 0;
}

static void decode_molecule(struct puzzle_molecule c, struct mechanism m, struct input_output *io)
{
    io->atoms = calloc(c.number_of_atoms, sizeof(io->atoms[0]));
    io->number_of_atoms = c.number_of_atoms;
    for (uint32_t i = 0; i < c.number_of_atoms; ++i) {
        io->atoms[i].atom = decode_atom(c.atoms[i].type);
        io->atoms[i].position = mechanism_relative_position(m, c.atoms[i].offset[0], c.atoms[i].offset[1], 1);
    }
    for (uint32_t i = 0; i < c.number_of_bonds; ++i) {
        struct puzzle_bond b = c.bonds[i];
        atom bonds = decode_bond_type(c.bonds[i].type);
        struct vector p1 = mechanism_relative_position(m, b.from[0], b.from[1], 1);
        struct vector p2 = mechanism_relative_position(m, b.to[0], b.to[1], 1);
        atom b1 = bonds & bond_direction(m, b.to[0] - b.from[0], b.to[1] - b.from[1]);
        atom b2 = bonds & bond_direction(m, b.from[0] - b.to[0], b.from[1] - b.to[1]);
        // yes, this is O(n^2).
        for (uint32_t j = 0; j < c.number_of_atoms; ++j) {
            struct vector p = io->atoms[j].position;
            if (p.u == p1.u && p.v == p1.v)
                io->atoms[j].atom |= b1;
            if (p.u == p2.u && p.v == p2.v)
                io->atoms[j].atom |= b2;
        }
    }
}

// thanks to Syx for this data
static struct vector cabinet_walls_Small[] = { {2,0}, {1,1}, {0,2}, {-1,2}, {-2,2}, {-2,1}, {-2,0}, {-1,-1}, {0,-2}, {1,-2}, {2,-2}, {2,-1} };
static struct vector cabinet_walls_SmallWide[] = { {3,0}, {2,1}, {1,2}, {0,2}, {-1,2}, {-2,2}, {-2,1}, {-2,0}, {-1,-1}, {0,-2}, {1,-2}, {2,-2}, {3,-2}, {3,-1} };
static struct vector cabinet_walls_SmallWider[] = { {4,0}, {3,1}, {2,2}, {1,2}, {0,2}, {-1,2}, {-2,2}, {-2,1}, {-2,0}, {-1,-1}, {0,-2}, {1,-2}, {2,-2}, {3,-2}, {4,-2}, {4,-1} };
static struct vector cabinet_walls_Medium[] = { {3,0}, {2,1}, {1,2}, {0,3}, {-1,3}, {-2,3}, {-3,3}, {-3,2}, {-3,1}, {-3,0}, {-2,-1}, {-1,-2}, {0,-3}, {1,-3}, {2,-3}, {3,-3}, {3,-2}, {3,-1} };
static struct vector cabinet_walls_MediumWide[] = { {4,0}, {3,1}, {2,2}, {1,3}, {0,3}, {-1,3}, {-2,3}, {-3,3}, {-3,2}, {-3,1}, {-3,0}, {-2,-1}, {-1,-2}, {0,-3}, {1,-3}, {2,-3}, {3,-3}, {4,-3}, {4,-2}, {4,-1} };
static struct vector cabinet_walls_Large[] = { {4,0}, {3,1}, {2,2}, {1,3}, {0,4}, {-1,4}, {-2,4}, {-3,4}, {-4,4}, {-4,3}, {-4,2}, {-4,1}, {-4,0}, {-3,-1}, {-2,-2}, {-1,-3}, {0,-4}, {1,-4}, {2,-4}, {3,-4}, {4,-4}, {4,-3}, {4,-2}, {4,-1} };

static struct vector cabinet_insides_Small[] = { {0,-1}, {1,-1}, {-1,0}, {0,0}, {1,0}, {-1,1}, {0,1} };
static struct vector cabinet_insides_SmallWide[] = { {0,-1}, {1,-1}, {2,-1}, {-1,0}, {0,0}, {1,0}, {2,0}, {-1,1}, {0,1}, {1,1} };
static struct vector cabinet_insides_SmallWider[] = { {0,-1}, {1,-1}, {2,-1}, {3,-1}, {-1,0}, {0,0}, {1,0}, {2,0}, {3,0}, {-1,1}, {0,1}, {1,1}, {2,1} };
static struct vector cabinet_insides_Medium[] = { {0,-2}, {1,-2}, {2,-2}, {-1,-1}, {0,-1}, {1,-1}, {2,-1}, {-2,0}, {-1,0}, {0,0}, {1,0}, {2,0}, {-2,1}, {-1,1}, {0,1}, {1,1}, {-2,2}, {-1,2}, {0,2} };
static struct vector cabinet_insides_MediumWide[] = { {0,-2}, {1,-2}, {2,-2}, {3,-2}, {-1,-1}, {0,-1}, {1,-1}, {2,-1}, {3,-1}, {-2,0}, {-1,0}, {0,0}, {1,0}, {2,0}, {3,0}, {-2,1}, {-1,1}, {0,1}, {1,1}, {2,1}, {-2,2}, {-1,2}, {0,2}, {1,2} };
static struct vector cabinet_insides_Large[] = { {0,-3}, {1,-3}, {2,-3}, {3,-3}, {-1,-2}, {0,-2}, {1,-2}, {2,-2}, {3,-2}, {-2,-1}, {-1,-1}, {0,-1}, {1,-1}, {2,-1}, {3,-1}, {-3,0}, {-2,0}, {-1,0}, {0,0}, {1,0}, {2,0}, {3,0}, {-3,1}, {-2,1}, {-1,1}, {0,1}, {1,1}, {2,1}, {-3,2}, {-2,2}, {-1,2}, {0,2}, {1,2}, {-3,3}, {-2,3}, {-1,3}, {0,3} };

static size_t number_of_walls_for_cabinet_type(struct byte_string type)
{
    if (byte_string_is(type, "Small"))
        return sizeof(cabinet_walls_Small) / sizeof(cabinet_walls_Small[0]);
    else if (byte_string_is(type, "SmallWide"))
        return sizeof(cabinet_walls_SmallWide) / sizeof(cabinet_walls_SmallWide[0]);
    else if (byte_string_is(type, "SmallWider"))
        return sizeof(cabinet_walls_SmallWider) / sizeof(cabinet_walls_SmallWider[0]);
    else if (byte_string_is(type, "Medium"))
        return sizeof(cabinet_walls_Medium) / sizeof(cabinet_walls_Medium[0]);
    else if (byte_string_is(type, "MediumWide"))
        return sizeof(cabinet_walls_MediumWide) / sizeof(cabinet_walls_MediumWide[0]);
    else if (byte_string_is(type, "Large"))
        return sizeof(cabinet_walls_Large) / sizeof(cabinet_walls_Large[0]);
    else
        return 0;
}

static size_t copy_walls_for_cabinet_type(struct byte_string type, struct vector *dest, int32_t u, int32_t v)
{
    struct vector *src;
    if (byte_string_is(type, "Small"))
        src = cabinet_walls_Small;
    else if (byte_string_is(type, "SmallWide"))
        src = cabinet_walls_SmallWide;
    else if (byte_string_is(type, "SmallWider"))
        src = cabinet_walls_SmallWider;
    else if (byte_string_is(type, "Medium"))
        src = cabinet_walls_Medium;
    else if (byte_string_is(type, "MediumWide"))
        src = cabinet_walls_MediumWide;
    else if (byte_string_is(type, "Large"))
        src = cabinet_walls_Large;
    else
        return 0;
    size_t n = number_of_walls_for_cabinet_type(type);
    for (size_t i = 0; i < n; ++i) {
        dest[i].u = src[i].u + u;
        dest[i].v = src[i].v + v;
    }
    return n;
}

static size_t number_of_insides_for_cabinet_type(struct byte_string type)
{
    if (byte_string_is(type, "Small"))
        return sizeof(cabinet_insides_Small) / sizeof(cabinet_insides_Small[0]);
    else if (byte_string_is(type, "SmallWide"))
        return sizeof(cabinet_insides_SmallWide) / sizeof(cabinet_insides_SmallWide[0]);
    else if (byte_string_is(type, "SmallWider"))
        return sizeof(cabinet_insides_SmallWider) / sizeof(cabinet_insides_SmallWider[0]);
    else if (byte_string_is(type, "Medium"))
        return sizeof(cabinet_insides_Medium) / sizeof(cabinet_insides_Medium[0]);
    else if (byte_string_is(type, "MediumWide"))
        return sizeof(cabinet_insides_MediumWide) / sizeof(cabinet_insides_MediumWide[0]);
    else if (byte_string_is(type, "Large"))
        return sizeof(cabinet_insides_Large) / sizeof(cabinet_insides_Large[0]);
    else
        return 0;
}

static struct vector* get_insides_for_cabinet_type(struct byte_string type)
{
    if (byte_string_is(type, "Small"))
        return cabinet_insides_Small;
    else if (byte_string_is(type, "SmallWide"))
        return cabinet_insides_SmallWide;
    else if (byte_string_is(type, "SmallWider"))
        return cabinet_insides_SmallWider;
    else if (byte_string_is(type, "Medium"))
        return cabinet_insides_Medium;
    else if (byte_string_is(type, "MediumWide"))
        return cabinet_insides_MediumWide;
    else if (byte_string_is(type, "Large"))
        return cabinet_insides_Large;
    else
        return 0;
}

static void mark_visible_region(struct solution *solution, struct vector p, int32_t hex_radius)
{
    if (p.u < INT32_MIN + hex_radius)
        solution->min_visible_u = INT32_MIN;
    else if (p.u - hex_radius < solution->min_visible_u)
        solution->min_visible_u = p.u - hex_radius;
    if (p.u > INT32_MAX - hex_radius)
        solution->max_visible_u = INT32_MAX;
    else if (p.u + hex_radius > solution->max_visible_u)
        solution->max_visible_u = p.u + hex_radius;
    if (p.v < INT32_MIN + hex_radius)
        solution->min_visible_v = INT32_MIN;
    else if (p.v - hex_radius < solution->min_visible_v)
        solution->min_visible_v = p.v - hex_radius;
    if (p.v > INT32_MAX - hex_radius)
        solution->max_visible_v = INT32_MAX;
    else if (p.v + hex_radius > solution->max_visible_v)
        solution->max_visible_v = p.v + hex_radius;
}

static void mark_visible_input_output(struct solution *solution, struct input_output *io)
{
    for (uint32_t i = 0; i < io->number_of_atoms; ++i)
        mark_visible_region(solution, io->atoms[i].position, 0);
}

#define MAP_SIZE 32

static uint8_t cabinet_for_position(uint8_t cabinet_map[MAP_SIZE][MAP_SIZE], struct vector position)
{
    int32_t u = position.u + (MAP_SIZE >> 1);
    int32_t v = position.v + (MAP_SIZE >> 1);
    if (u >= 0 && u < MAP_SIZE && v >= 0 && v < MAP_SIZE)
        return cabinet_map[u][v];
    return 0;
}

static struct vector decode_position(signed char position[2])
{
    return (struct vector) { position[0], position[1] };
}

static char* check_production_constraints(struct solution *solution, struct puzzle_production_info *info)
{
    // map u,v coordinates => cabinet id + 1 (0 means no cabinet)
    uint8_t cabinet_map[MAP_SIZE][MAP_SIZE] = { 0 };
    for (uint32_t i = 0; i < info->number_of_cabinets; ++i) {
        struct vector *src = get_insides_for_cabinet_type(info->cabinets[i].type);
        size_t n = number_of_insides_for_cabinet_type(info->cabinets[i].type);
        for (size_t j = 0; j < n; ++j) {
            int32_t u = info->cabinets[i].position[0] + src[j].u + (MAP_SIZE >> 1);
            int32_t v = info->cabinets[i].position[1] + src[j].v + (MAP_SIZE >> 1);
            if (u >= 0 && u < MAP_SIZE && v >= 0 && v < MAP_SIZE)
                cabinet_map[u][v] = i + 1;
        }
    }

    // check conduits
    bool swapped = false;
    for (uint32_t i = 0; i < solution->number_of_conduits; ++i) {
        struct conduit *conduit = &solution->conduits[i];
        uint32_t conduit_index = conduit->id - 100; // conduit ids start at 100
        if (conduit_index >= info->number_of_conduits)
            return "solution contains a conduit not defined in the puzzle file";

        struct puzzle_conduit *puzzle_conduit = &info->conduits[conduit_index];
        if (conduit->number_of_positions != puzzle_conduit->number_of_hexes)
            return "solution contains a conduit with an edited shape";

        signed char *starting_position = (i & 1) ? puzzle_conduit->starting_position_b : puzzle_conduit->starting_position_a;
        uint8_t expected_cabinet = cabinet_for_position(cabinet_map, decode_position(starting_position));

        // determine if this conduit is the A-side or the B-side of its pair
        if ((i & 1) == 0) {
            struct vector center = solution->glyphs[conduit->glyph_index].position;
            swapped = cabinet_for_position(cabinet_map, center) != expected_cabinet;
        }
        struct mechanism *m = &solution->glyphs[swapped ? conduit->other_side_glyph_index : conduit->glyph_index];

        for (uint32_t j = 0; j < conduit->number_of_positions; ++j) {
            struct vector pos = conduit->positions[j];
            if (!vectors_equal(pos, decode_position(puzzle_conduit->hexes[j].offset)))
                return "solution contains a conduit with an edited shape";
            pos = mechanism_relative_position(*m, pos.u, pos.v, 1);
            if (cabinet_for_position(cabinet_map, pos) != expected_cabinet)
                return "solution moved a conduit outside of its original chamber";
        }
    }

    // check glyphs
    for (uint32_t i = 0; i < solution->number_of_glyphs; ++i) {
        struct mechanism *m = &solution->glyphs[i];
        if (m->type == CONDUIT)
            continue; // conduits were already checked
        const struct vector *footprint = glyph_footprint(m->type);
        for (int j = 0; ; ++j) {
            struct vector pos = mechanism_relative_position(*m, footprint[j].u, footprint[j].v, 1);
            if (cabinet_for_position(cabinet_map, pos) == 0)
                return "solution contains a glyph outside of the cabinet";
            if (vectors_equal(footprint[j], zero_vector))
                break;
        }
    }

    // check arms
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        struct mechanism *m = &solution->arms[i];
        uint8_t base_cabinet = cabinet_for_position(cabinet_map, m->position);
        if (base_cabinet == 0)
            return "solution contains an arm outside of the cabinet";
        int step = angular_distance_between_grabbers(m->type);
        for (int direction = 0; direction < 6; direction += step) {
            struct vector offset = u_offset_for_direction(direction);
            struct vector pos = mechanism_relative_position(*m, offset.u, offset.v, 1);
            uint8_t cabinet = cabinet_for_position(cabinet_map, pos);
            if (cabinet == 0)
                return "solution contains an arm outside of the cabinet";
            else if (cabinet != base_cabinet)
                return "solution contains an arm reaching across cabinet walls";
        }
    }

    // check tracks
    for (uint32_t i = 0; i < solution->track_table_size; ++i) {
        struct vector position = solution->track_positions[i];
        if (position.u == INT32_MIN && position.v == INT32_MIN)
            continue;
        if (cabinet_for_position(cabinet_map, position) == 0)
            return "solution contains a track outside of the cabinet";
    }

    // check inputs/outputs
    for (uint32_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
        struct input_output *io = &solution->inputs_and_outputs[i];
        for (uint32_t j = 0; j < io->number_of_atoms; ++j)
            if (cabinet_for_position(cabinet_map, io->atoms[j].position) == 0)
                return "solution contains an input/output outside of the cabinet";
    }

    // check isolation
    if (info->isolate_inputs_from_outputs) {
        enum input_output_type *isolation_data = calloc(info->number_of_cabinets + 1, sizeof(enum input_output_type));
        for (uint32_t i = 0; i < solution->number_of_inputs_and_outputs; ++i) {
            struct input_output *io = &solution->inputs_and_outputs[i];
            uint8_t cabinet = cabinet_for_position(cabinet_map, io->atoms->position);
            isolation_data[cabinet] |= io->type;
            if ((isolation_data[cabinet] & INPUT) && (isolation_data[cabinet] & OUTPUT))
                return "solution breaks the input/output isolation constraint";
        }
        free(isolation_data);
    }

    return NULL;
}

bool decode_solution(struct solution *solution, struct puzzle_file *pf, struct solution_file *sf, const char **error)
{
    const char *ignored_error;
    if (!error)
        error = &ignored_error;
    for (uint32_t i = 0; i < pf->number_of_inputs; ++i) {
        if (pf->inputs[i].number_of_atoms == 0) {
            *error = "puzzle file contains a reagent with no atoms";
            return false;
        }
    }
    int64_t number_of_unplaced_outputs = 0;
    for (uint32_t i = 0; i < pf->number_of_outputs; ++i) {
        number_of_unplaced_outputs++;
        if (pf->outputs[i].number_of_atoms == 0) {
            *error = "puzzle file contains a product with no atoms";
            return false;
        }
    }

    size_t number_of_arms = 0;
    size_t number_of_glyphs = 0;
    size_t number_of_conduits = 0;
    size_t number_of_track_hexes = 0;
    size_t number_of_inputs_and_outputs = 0;
    // first pass through the solution file: count how many things of each type
    // there are.  these counts are used to allocate arrays of the correct size.
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        uint64_t parts_available = parts_available_bits_for_part_name(sf->parts[i].name);
        if ((parts_available & pf->parts_available) != parts_available) {
            *error = "solution contains a part that has been disabled in the puzzle file";
            return false;
        }
        enum mechanism_type type = decode_mechanism_type(sf->parts[i].name);
        if (type & ANY_ARM) {
            if (sf->parts[i].size > 3) {
                *error = "solution contains a too-long arm";
                return false;
            } else if (sf->parts[i].size < 1) {
                *error = "solution contains a zero-length arm";
                return false;
            }
            number_of_arms++;
        } else if (type & ANY_GLYPH)
            number_of_glyphs++;
        else if (byte_string_is(sf->parts[i].name, "track"))
            number_of_track_hexes += sf->parts[i].number_of_track_hexes;
        else if (byte_string_is(sf->parts[i].name, "input")) {
            if (sf->parts[i].which_input_or_output >= pf->number_of_inputs) {
                *error = "solution refers to an input that doesn't exist in the puzzle file";
                return false;
            }
            number_of_inputs_and_outputs++;
        } else if (byte_string_is(sf->parts[i].name, "out-std")) {
            if (sf->parts[i].which_input_or_output >= pf->number_of_outputs) {
                *error = "solution refers to an output that doesn't exist in the puzzle file";
                return false;
            }
            number_of_inputs_and_outputs++;
            number_of_unplaced_outputs--;
        } else if (byte_string_is(sf->parts[i].name, "out-rep")) {
            uint32_t which_output = sf->parts[i].which_input_or_output;
            if (which_output >= pf->number_of_outputs) {
                *error = "solution refers to an output that doesn't exist in the puzzle file";
                return false;
            }
            number_of_inputs_and_outputs++;
            number_of_unplaced_outputs--;
        } else if (byte_string_is(sf->parts[i].name, "pipe")) {
            number_of_glyphs++;
            number_of_conduits++;
        }
    }
    if (number_of_unplaced_outputs != 0) {
        *error = "all products in the puzzle must be placed in the solution";
        destroy(solution, 0);
        return false;
    }
    solution->number_of_arms = number_of_arms;
    solution->number_of_glyphs = number_of_glyphs;
    solution->number_of_conduits = number_of_conduits;
    solution->number_of_inputs_and_outputs = number_of_inputs_and_outputs;

    size_t number_of_cabinet_walls = 0;
    for (uint32_t i = 0; pf->production_info && i < pf->production_info->number_of_cabinets; ++i)
        number_of_cabinet_walls += number_of_walls_for_cabinet_type(pf->production_info->cabinets[i].type);
    solution->number_of_cabinet_walls = number_of_cabinet_walls;
    solution->cabinet_walls = calloc(number_of_cabinet_walls, sizeof(struct vector));
    number_of_cabinet_walls = 0;
    for (uint32_t i = 0; pf->production_info && i < pf->production_info->number_of_cabinets; ++i) {
        number_of_cabinet_walls += copy_walls_for_cabinet_type(pf->production_info->cabinets[i].type,
            solution->cabinet_walls + number_of_cabinet_walls,
            pf->production_info->cabinets[i].position[0],
            pf->production_info->cabinets[i].position[1]);
    }

    // now that we know how many elements each array should have, allocate them
    // all here.
    solution->glyphs = calloc(solution->number_of_glyphs, sizeof(struct mechanism));

    solution->arms = calloc(solution->number_of_arms, sizeof(struct mechanism));
    solution->arm_tape = calloc(solution->number_of_arms, sizeof(char *));
    solution->arm_tape_length = calloc(solution->number_of_arms, sizeof(size_t));
    solution->arm_tape_start_cycle = calloc(solution->number_of_arms, sizeof(int64_t));

    solution->conduits = calloc(solution->number_of_conduits, sizeof(struct conduit));

    solution->inputs_and_outputs = calloc(solution->number_of_inputs_and_outputs, sizeof(struct input_output));

    solution->track_table_size = 1;
    while (2 * solution->track_table_size < 3 * number_of_track_hexes)
        solution->track_table_size *= 2;
    solution->track_positions = calloc(solution->track_table_size, sizeof(solution->track_positions[0]));
    for (int i = 0; i < solution->track_table_size; ++i)
        solution->track_positions[i] = (struct vector){ INT32_MIN, INT32_MIN };
    solution->track_plus_motions = calloc(solution->track_table_size, sizeof(solution->track_plus_motions[0]));
    solution->track_minus_motions = calloc(solution->track_table_size, sizeof(solution->track_minus_motions[0]));

    solution->min_visible_u = INT32_MAX;
    solution->max_visible_u = INT32_MIN;
    solution->min_visible_v = INT32_MAX;
    solution->max_visible_v = INT32_MIN;

    // second pass: fill in the arrays with the data from the file.  this pass
    // goes in reverse to properly handle overlapping tracks.
    uint32_t arm_index = solution->number_of_arms - 1;
    uint32_t glyph_index = solution->number_of_glyphs - 1;
    uint32_t conduit_index = solution->number_of_conduits - 1;
    size_t io_index = solution->number_of_inputs_and_outputs - 1;
    for (uint32_t i = sf->number_of_parts - 1; i < sf->number_of_parts; --i) {
        struct solution_part part = sf->parts[i];
        part.rotation = part.rotation % 6;
        struct mechanism m = {
            .type = decode_mechanism_type(part.name),
            .position = { part.position[0], part.position[1] },
            .direction_u = u_offset_for_direction(part.rotation),
            .direction_v = v_offset_for_direction(part.rotation),
            .arm_rotation = part.rotation,
        };
        if (m.type & ANY_ARM) {
            m.direction_u.u *= part.size;
            m.direction_u.v *= part.size;
            m.direction_v.u *= part.size;
            m.direction_v.v *= part.size;
            solution->arms[arm_index--] = m;
            mark_visible_region(solution, m.position, (m.type & PISTON) ? 3 : part.size);
        } else if (m.type & ANY_GLYPH) {
            solution->glyphs[glyph_index--] = m;
            if (!(m.type & EQUILIBRIUM)) {
                const struct vector *footprint = glyph_footprint(m.type);
                for (int j = 0; ; j++) {
                    struct vector p = footprint[j];
                    mark_visible_region(solution, mechanism_relative_position(m, p.u, p.v, 1), 0);
                    if (vectors_equal(p, zero_vector))
                        break;
                }
            }
        } else if (byte_string_is(part.name, "track")) {
            struct vector last_position = m.position;
            for (uint32_t j = 0; part.number_of_track_hexes > 0 && j < part.number_of_track_hexes + 1; ++j) {
                struct solution_hex_offset hex;
                if (j < part.number_of_track_hexes)
                    hex = part.track_hexes[j];
                else {
                    hex = part.track_hexes[0];
                    int32_t du = hex.offset[0] - part.track_hexes[j - 1].offset[0];
                    int32_t dv = hex.offset[1] - part.track_hexes[j - 1].offset[1];
                    // two-hex tracks can't become a loop.  also, if the offset
                    // between the two hexes isn't a cardinal direction, then
                    // the ends are too far away for the track to become a loop.
                    if (part.number_of_track_hexes <= 2 || direction_for_offset((struct vector){ du, dv }) < 0) {
                        // ensure a plus motion leaves the arm in place.
                        uint32_t index;
                        lookup_track(solution, last_position, &index);
                        solution->track_plus_motions[index] = (struct vector){ 0, 0 };
                        break;
                    }
                }
                struct vector p = mechanism_relative_position(m, hex.offset[0], hex.offset[1], 1);
                if (j == 0)
                    last_position = p;
                uint32_t index;
                lookup_track(solution, p, &index);
                if (j < part.number_of_track_hexes && (solution->track_positions[index].u != INT32_MIN || solution->track_positions[index].v != INT32_MIN))
                    solution->track_self_overlap++;
                solution->track_positions[index] = p;
                solution->track_minus_motions[index] = (struct vector){ last_position.u - p.u, last_position.v - p.v };
                lookup_track(solution, last_position, &index);
                solution->track_plus_motions[index] = (struct vector){ p.u - last_position.u, p.v - last_position.v };
                last_position = p;
                mark_visible_region(solution, p, 3);
            }
        } else if (byte_string_is(part.name, "pipe")) {
            m.type = CONDUIT;
            struct conduit conduit = {
                .glyph_index = glyph_index,
                .id = part.conduit_id,
                .number_of_positions = part.number_of_conduit_hexes,
                .positions = calloc(sizeof(struct vector), part.number_of_conduit_hexes),
                .atoms = calloc(sizeof(struct atom_at_position), part.number_of_conduit_hexes),
                .molecule_lengths = calloc(sizeof(uint32_t), part.number_of_conduit_hexes),
            };
            for (uint32_t j = 0; j < part.number_of_conduit_hexes; ++j) {
                conduit.positions[j].u = part.conduit_hexes[j].offset[0];
                conduit.positions[j].v = part.conduit_hexes[j].offset[1];
                mark_visible_region(solution, mechanism_relative_position(m, conduit.positions[j].u, conduit.positions[j].v, 1), 0);
            }
            solution->glyphs[glyph_index--] = m;
            solution->conduits[conduit_index--] = conduit;
        } else if (byte_string_is(part.name, "input")) {
            struct puzzle_molecule c = pf->inputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = INPUT;
            io->puzzle_index = part.which_input_or_output;
            io->solution_index = i;
            decode_molecule(c, m, io);
            io_index--;
            mark_visible_input_output(solution, io);
        } else if (byte_string_is(part.name, "out-std")) {
            struct puzzle_molecule c = pf->outputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = SINGLE_OUTPUT;
            io->puzzle_index = part.which_input_or_output;
            io->solution_index = i;
            decode_molecule(c, m, io);
            io_index--;
            mark_visible_input_output(solution, io);
        } else if (byte_string_is(part.name, "out-rep")) {
            struct puzzle_molecule c = pf->outputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = REPEATING_OUTPUT;
            io->puzzle_index = part.which_input_or_output;
            io->solution_index = i;
            decode_molecule(c, m, io);
            io->original_atoms = io->atoms;
            io->number_of_original_atoms = io->number_of_atoms;
            if (io->number_of_original_atoms == 0) {
                *error = "solution contains an empty infinite product";
                destroy(solution, 0);
                return false;
            }
            struct atom_at_position *placeholder = &io->original_atoms[io->number_of_original_atoms - 1];
            for (uint32_t i = 0; i < io->number_of_original_atoms; ++i) {
                if (!(io->original_atoms[i].atom & REPEATING_OUTPUT_PLACEHOLDER))
                    continue;
                struct atom_at_position a = io->original_atoms[i];
                io->original_atoms[i] = *placeholder;
                *placeholder = a;
                break;
            }
            if (!(placeholder->atom & REPEATING_OUTPUT_PLACEHOLDER)) {
                *error = "solution contains an infinite product without a repetition placeholder";
                destroy(solution, 0);
                return false;
            }
            io->atoms = 0;
            io->number_of_atoms = 0;
            io->repetition_origin = m.position;
            io->outputs_per_repetition = pf->output_scale;
            if (!repeat_molecule(io, REPEATING_OUTPUT_REPETITIONS, error)) {
                destroy(solution, 0);
                return false;
            }
            io_index--;
            mark_visible_input_output(solution, io);
        }
    }
    // sort conduits by id in order to find pairs of linked conduits.
    qsort(solution->conduits, solution->number_of_conduits,
     sizeof(struct conduit), compare_conduits_by_id);
    struct conduit *destination = 0;
    for (uint32_t i = 0; i < solution->number_of_conduits; ++i) {
        struct conduit *conduit = &solution->conduits[i];
        struct conduit *next = 0;
        if (i + 1 < solution->number_of_conduits)
            next = &solution->conduits[i + 1];
        if (destination && destination->id == conduit->id)
            conduit->other_side_glyph_index = destination->glyph_index;
        else if (!next || next->id != conduit->id) {
            *error = "solution contains an unpaired conduit";
            destroy(solution, 0);
            return false;
        } else {
            destination = conduit;
            conduit->other_side_glyph_index = next->glyph_index;
        }
        solution->glyphs[conduit->glyph_index].conduit_index = i;
    }
    // production-only pass: enforce production constraints
    if (pf->production_info) {
        *error = check_production_constraints(solution, pf->production_info);
        if (*error)
            return false;
    } else if (solution->number_of_conduits > 0) {
        *error = "solution contains a conduit not defined in the puzzle file";
        return false;
    }
    // decode arm tapes in one final pass.  this has to be another pass because
    // reset instructions depend on where track has been placed.
    arm_index = 0;
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        struct solution_part part = sf->parts[i];
        if (!(decode_mechanism_type(part.name) & ANY_ARM))
            continue;
        if (!part.number_of_instructions) {
            arm_index++;
            continue;
        }
        qsort(part.instructions, part.number_of_instructions,
         sizeof(part.instructions[0]), compare_instructions_by_index);
        int32_t min_tape = part.instructions[0].index;
        int32_t max_tape = part.instructions[part.number_of_instructions - 1].index;
        if ((int64_t)max_tape - (int64_t)min_tape > 99999) {
            *error = "solution has an arm with an instruction tape that's too long";
            destroy(solution, 0);
            return false;
        }
        int32_t tape_length = max_tape - min_tape + 1;
        if ((int64_t)min_tape + (int64_t)tape_length * 2 > INT32_MAX) {
            *error = "solution has an arm with a potential instruction tape index overflow";
            destroy(solution, 0);
            return false;
        }
        // multiply by two to leave room for reset and repeat instructions
        // (which cannot exceed the length of the earlier instructions even in
        // the worst case).
        solution->arm_tape[arm_index] = calloc(tape_length * 2, 1);
        solution->arm_tape_start_cycle[arm_index] = min_tape;
        int32_t last_end = 0;
        int32_t last_repeat = 0;
        int32_t reset_from = 0;
        char *tape = solution->arm_tape[arm_index];
        for (uint32_t j = 0; j < part.number_of_instructions; ++j) {
            struct solution_instruction inst = part.instructions[j];
            if (j > 0 && inst.index == part.instructions[j - 1].index) {
                *error = "solution contains an arm with two instructions that have the same index";
                destroy(solution, 0);
                return false;
            }
            uint64_t parts_available = parts_available_bits_for_instruction(inst.instruction);
            if ((parts_available & pf->parts_available) != parts_available) {
                *error = "solution contains an instruction that has been disabled in the puzzle file";
                destroy(solution, 0);
                return false;
            }
            int32_t n = inst.index - min_tape;
            if (inst.instruction == 'C') { // repeat
                if (last_repeat < -min_tape)
                    last_repeat = -min_tape;
                while (j < part.number_of_instructions && part.instructions[j].instruction == 'C') {
                    if (last_end > part.instructions[j].index - min_tape) {
                        *error = "solution contains a repeat instruction that overlaps with a reset instruction";
                        destroy(solution, 0);
                        return false;
                    }
                    if (last_end > last_repeat) {
                        memcpy(tape + part.instructions[j].index - min_tape, tape + last_repeat, last_end - last_repeat);
                        int32_t m = part.instructions[j].index - min_tape + last_end - last_repeat;
                        if (m > tape_length)
                            tape_length = m;
                    }
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
                int track_looping_steps = 0;
                int rotation = 0;
                int piston = part.size;
                int grab = 0;
                for (uint32_t k = reset_from; k < n; ++k) {
                    struct vector motion = { 0, 0 };
                    if (tape[k] == 'a') {
                        rotation++;
                    } else if (tape[k] == 'd') {
                        rotation--;
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
                    if (track.u == position.u && track.v == position.v) {
                        track_looping_steps += track_steps;
                        track_steps = 0;
                    }
                }
                if (grab > 0)
                    tape[n++] = 'f';
                while (piston > part.size) {
                    tape[n++] = 's';
                    piston--;
                }
                while (rotation > 3)
                    rotation -= 6;
                while (rotation < -3)
                    rotation += 6;
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
                    // if taking into account looping steps changes the sign,
                    // resolve ties in the opposite way.
                    if (track_steps * (track_steps + track_looping_steps) < 0)
                        search_depth++;
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
                while (piston < part.size) {
                    tape[n++] = 'w';
                    piston++;
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
    int64_t tape_start_cycle = INT64_MAX;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        if (solution->arm_tape_start_cycle[i] < tape_start_cycle)
            tape_start_cycle = solution->arm_tape_start_cycle[i];
    }
    for (uint32_t i = 0; i < solution->number_of_arms; ++i)
        solution->arm_tape_start_cycle[i] -= tape_start_cycle;

    solution->target_number_of_outputs = 6 * pf->output_scale;
    return true;
}

uint64_t solution_file_cost(struct solution_file *sf)
{
    uint64_t cost = 0;
    for (uint32_t i = 0; i < sf->number_of_parts; ++i) {
        struct byte_string part_name = sf->parts[i].name;
        if (byte_string_is(part_name, "glyph-calcification"))
            cost += 10;
        else if (byte_string_is(part_name, "glyph-life-and-death"))
            cost += 20;
        else if (byte_string_is(part_name, "glyph-projection"))
            cost += 20;
        else if (byte_string_is(part_name, "glyph-dispersion"))
            cost += 20;
        else if (byte_string_is(part_name, "glyph-purification"))
            cost += 20;
        else if (byte_string_is(part_name, "glyph-duplication"))
            cost += 20;
        else if (byte_string_is(part_name, "glyph-unification"))
            cost += 20;
        else if (byte_string_is(part_name, "bonder"))
            cost += 10;
        else if (byte_string_is(part_name, "unbonder"))
            cost += 10;
        else if (byte_string_is(part_name, "bonder-prisma"))
            cost += 20;
        else if (byte_string_is(part_name, "bonder-speed"))
            cost += 30;
        else if (byte_string_is(part_name, "arm1"))
            cost += 20;
        else if (byte_string_is(part_name, "arm2"))
            cost += 30;
        else if (byte_string_is(part_name, "arm3"))
            cost += 30;
        else if (byte_string_is(part_name, "arm6"))
            cost += 30;
        else if (byte_string_is(part_name, "piston"))
            cost += 40;
        else if (byte_string_is(part_name, "baron"))
            cost += 30;
        else if (byte_string_is(part_name, "track"))
            cost += sf->parts[i].number_of_track_hexes * 5;
    }
    return cost;
}

uint64_t solution_instructions(struct solution *solution)
{
    uint64_t instructions = 0;
    for (uint32_t i = 0; i < solution->number_of_arms; ++i) {
        for (uint32_t j = 0; j < solution->arm_tape_length[i]; ++j) {
            if (solution->arm_tape[i][j] != ' ' && solution->arm_tape[i][j] != '\0')
                instructions++;
        }
    }
    return instructions;
}
