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
        io->atoms[i].position = mechanism_relative_position(m, c.atoms[i].offset[1], c.atoms[i].offset[0], 1);
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
            struct vector p = io->atoms[j].position;
            if (p.u == p1.u && p.v == p1.v)
                io->atoms[j].atom |= b1;
            if (p.u == p2.u && p.v == p2.v)
                io->atoms[j].atom |= b2;
        }
    }
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

    // second pass: fill in the arrays with the data from the file.  this pass
    // goes in reverse to properly handle overlapping tracks.
    uint32_t arm_index = solution->number_of_arms - 1;
    uint32_t glyph_index = solution->number_of_glyphs - 1;
    uint32_t conduit_index = solution->number_of_conduits - 1;
    size_t io_index = solution->number_of_inputs_and_outputs - 1;
    for (uint32_t i = sf->number_of_parts - 1; i < sf->number_of_parts; --i) {
        struct solution_part part = sf->parts[i];
        struct mechanism m = {
            .type = decode_mechanism_type(part.name),
            .position = { part.position[1], part.position[0] },
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
        } else if (m.type & ANY_GLYPH)
            solution->glyphs[glyph_index--] = m;
        else if (byte_string_is(part.name, "track")) {
            struct vector last_position = m.position;
            for (uint32_t j = 0; j < part.number_of_track_hexes + 1; ++j) {
                struct solution_hex_offset hex;
                if (j < part.number_of_track_hexes)
                    hex = part.track_hexes[j];
                else {
                    hex = part.track_hexes[0];
                    int32_t du = hex.offset[1] - part.track_hexes[j - 1].offset[1];
                    int32_t dv = hex.offset[0] - part.track_hexes[j - 1].offset[0];
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
                struct vector p = mechanism_relative_position(m, hex.offset[1], hex.offset[0], 1);
                if (j == 0)
                    last_position = p;
                uint32_t index;
                lookup_track(solution, p, &index);
                solution->track_positions[index] = p;
                solution->track_minus_motions[index] = (struct vector){ last_position.u - p.u, last_position.v - p.v };
                lookup_track(solution, last_position, &index);
                solution->track_plus_motions[index] = (struct vector){ p.u - last_position.u, p.v - last_position.v };
                last_position = p;
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
                conduit.positions[j].u = part.conduit_hexes[j].offset[1];
                conduit.positions[j].v = part.conduit_hexes[j].offset[0];
            }
            solution->glyphs[glyph_index--] = m;
            solution->conduits[conduit_index--] = conduit;
        } else if (byte_string_is(part.name, "input")) {
            struct puzzle_molecule c = pf->inputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = INPUT;
            io->puzzle_index = part.which_input_or_output;
            decode_molecule(c, m, io);
            io_index--;
        } else if (byte_string_is(part.name, "out-std")) {
            struct puzzle_molecule c = pf->outputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = SINGLE_OUTPUT;
            io->puzzle_index = part.which_input_or_output;
            decode_molecule(c, m, io);
            io_index--;
        } else if (byte_string_is(part.name, "out-rep")) {
            struct puzzle_molecule c = pf->outputs[part.which_input_or_output];
            struct input_output *io = &solution->inputs_and_outputs[io_index];
            io->type = REPEATING_OUTPUT;
            io->puzzle_index = part.which_input_or_output;
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
            int32_t n = inst.index - min_tape;
            if (inst.instruction == 'C') { // repeat
                while (j < part.number_of_instructions && part.instructions[j].instruction == 'C') {
                    if (last_end > part.instructions[j].index - min_tape) {
                        *error = "solution contains a repeat instruction that overlaps with a reset instruction";
                        destroy(solution, 0);
                        return false;
                    }
                    memcpy(tape + part.instructions[j].index - min_tape, tape + last_repeat, last_end - last_repeat);
                    int32_t m = part.instructions[j].index - min_tape + last_end - last_repeat;
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
                    if (track.u == position.u && track.v == position.v)
                        track_steps = 0;
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
