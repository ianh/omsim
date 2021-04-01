#include "parse.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

static struct byte_string read_file(const char *filename)
{
    struct byte_string contents;
    FILE *f = fopen(filename, "r");
    fseek(f, 0, SEEK_END);
    contents.length = ftell(f);
    contents.bytes = malloc(contents.length);
    fseek(f, 0, SEEK_SET);
    fread(contents.bytes, contents.length, 1, f);
    fclose(f);
    return contents;
}

static unsigned char read_byte(struct byte_string *contents)
{
    if (contents->length < 1)
        return 0;
    unsigned char value = contents->bytes[0];
    contents->length--;
    contents->bytes++;
    return value;
}

static signed char read_signed_byte(struct byte_string *contents)
{
    // this is probably overly pedantic, but it avoids signed overflow.
    unsigned char value = read_byte(contents);
    if (value <= SCHAR_MAX)
        return value;
    return -(signed char)((unsigned char)~value)-1;
}

static uint32_t read_uint32(struct byte_string *contents)
{
    if (contents->length < 4)
        return 0;
    uint32_t value = contents->bytes[0] | (contents->bytes[1] << 8) |
     (contents->bytes[2] << 16) | (contents->bytes[3] << 24);
    contents->length -= 4;
    contents->bytes += 4;
    return value;
}

static int32_t read_int32(struct byte_string *contents)
{
    uint32_t value = read_uint32(contents);
    if (value <= INT32_MAX)
        return value;
    return -(int32_t)((uint32_t)~value)-1;
}

static int32_t read_uint64(struct byte_string *contents)
{
    uint32_t lo = read_uint32(contents);
    uint32_t hi = read_uint32(contents);
    return (uint64_t)lo | ((uint64_t)hi << 32);
}

static struct byte_string read_string(struct byte_string *contents)
{
    size_t len = 0;
    int shift = 0;
    while (contents->length > 0) {
        unsigned char byte = read_byte(contents);
        len |= (byte & 0x7f) << shift;
        shift += 7;
        if (!(byte & 0x80))
            break;
    }
    if (contents->length < len)
        len = contents->length;
    struct byte_string value = { contents->bytes, len };
    contents->length -= len;
    contents->bytes += len;
    return value;
}

static void parse_puzzle_molecule(struct byte_string *b, struct puzzle_molecule *molecule)
{
    molecule->number_of_atoms = read_uint32(b);
    molecule->atoms = calloc(molecule->number_of_atoms, sizeof(struct puzzle_atom));
    for (uint32_t j = 0; j < molecule->number_of_atoms; ++j) {
        molecule->atoms[j].type = read_byte(b);
        molecule->atoms[j].offset[0] = read_signed_byte(b);
        molecule->atoms[j].offset[1] = read_signed_byte(b);
    }
    molecule->number_of_bonds = read_uint32(b);
    molecule->bonds = calloc(molecule->number_of_bonds, sizeof(struct puzzle_bond));
    for (uint32_t j = 0; j < molecule->number_of_bonds; ++j) {
        molecule->bonds[j].type = read_byte(b);
        molecule->bonds[j].from[0] = read_signed_byte(b);
        molecule->bonds[j].from[1] = read_signed_byte(b);
        molecule->bonds[j].to[0] = read_signed_byte(b);
        molecule->bonds[j].to[1] = read_signed_byte(b);
    }
}

struct puzzle_file *parse_puzzle_byte_string(struct byte_string b)
{
    struct puzzle_file *puzzle = calloc(1, sizeof(struct puzzle_file));
    puzzle->bytes = b.bytes;
    if (read_uint32(&b) != 3) {
        free_puzzle_file(puzzle);
        return 0;
    }
    puzzle->name = read_string(&b);
    puzzle->creator = read_uint64(&b);
    puzzle->parts_available = read_uint64(&b);
    puzzle->number_of_inputs = read_uint32(&b);
    puzzle->inputs = calloc(puzzle->number_of_inputs, sizeof(struct puzzle_molecule));
    for (uint32_t i = 0; i < puzzle->number_of_inputs; ++i)
        parse_puzzle_molecule(&b, &puzzle->inputs[i]);
    puzzle->number_of_outputs = read_uint32(&b);
    puzzle->outputs = calloc(puzzle->number_of_outputs, sizeof(struct puzzle_molecule));
    for (uint32_t i = 0; i < puzzle->number_of_outputs; ++i)
        parse_puzzle_molecule(&b, &puzzle->outputs[i]);
    puzzle->output_scale = read_uint32(&b);
    // xx support production puzzles
    return puzzle;
}

struct puzzle_file *parse_puzzle_file(const char *path)
{
    struct byte_string b = read_file(path);
    struct puzzle_file *puzzle = parse_puzzle_byte_string(b);
    if (!puzzle)
        free(b.bytes);
    puzzle->owns_bytes = true;
    return puzzle;
}

void free_puzzle_file(struct puzzle_file *puzzle)
{
    for (uint32_t i = 0; i < puzzle->number_of_inputs; ++i) {
        free(puzzle->inputs[i].atoms);
        free(puzzle->inputs[i].bonds);
    }
    free(puzzle->inputs);
    for (uint32_t i = 0; i < puzzle->number_of_outputs; ++i) {
        free(puzzle->outputs[i].atoms);
        free(puzzle->outputs[i].bonds);
    }
    free(puzzle->outputs);
    if (puzzle->owns_bytes)
        free(puzzle->bytes);
    free(puzzle);
}

struct solution_file *parse_solution_byte_string(struct byte_string b)
{
    struct solution_file *solution = calloc(1, sizeof(struct solution_file));
    solution->bytes = b.bytes;
    if (read_uint32(&b) != 7) {
        free_solution_file(solution);
        return 0;
    }
    solution->puzzle = read_string(&b);
    solution->name = read_string(&b);
    if (read_uint32(&b)) {
        uint32_t zero = read_uint32(&b);
        solution->cycles = read_uint32(&b);
        uint32_t one = read_uint32(&b);
        solution->cost = read_uint32(&b);
        uint32_t two = read_uint32(&b);
        solution->area = read_uint32(&b);
        uint32_t three = read_uint32(&b);
        solution->instructions = read_uint32(&b);
        if (zero != 0 || one != 1 || two != 2 || three != 3) {
            free_solution_file(solution);
            return 0;
        }
    }
    solution->number_of_parts = read_uint32(&b);
    solution->parts = calloc(solution->number_of_parts, sizeof(struct solution_part));
    for (uint32_t i = 0; i < solution->number_of_parts; ++i) {
        struct solution_part *part = &solution->parts[i];
        part->name = read_string(&b);
        if (read_byte(&b) != 1) {
            free_solution_file(solution);
            return 0;
        }
        part->position[0] = read_int32(&b);
        part->position[1] = read_int32(&b);
        part->size = read_int32(&b);
        part->rotation = read_int32(&b);
        part->which_input_or_output = read_uint32(&b);

        part->number_of_instructions = read_uint32(&b);
        part->instructions = calloc(part->number_of_instructions, sizeof(struct solution_instruction));
        for (uint32_t j = 0; j < part->number_of_instructions; ++j) {
            part->instructions[j].index = read_int32(&b);
            part->instructions[j].instruction = read_byte(&b);
        }
        if (byte_string_is(part->name, "track")) {
            part->number_of_track_hexes = read_uint32(&b);
            part->track_hexes = calloc(part->number_of_track_hexes, sizeof(struct solution_hex_offset));
            for (uint32_t j = 0; j < part->number_of_track_hexes; ++j) {
                part->track_hexes[j].offset[0] = read_int32(&b);
                part->track_hexes[j].offset[1] = read_int32(&b);
            }
        }
        part->arm_number = read_uint32(&b);
        if (byte_string_is(part->name, "pipe")) {
            part->conduit_id = read_uint32(&b);
            part->number_of_conduit_hexes = read_uint32(&b);
            part->conduit_hexes = calloc(part->number_of_conduit_hexes, sizeof(struct solution_hex_offset));
            for (uint32_t j = 0; j < part->number_of_conduit_hexes; ++j) {
                part->conduit_hexes[j].offset[0] = read_int32(&b);
                part->conduit_hexes[j].offset[1] = read_int32(&b);
            }
        }
    }
    return solution;
}

struct solution_file *parse_solution_file(const char *path)
{
    struct byte_string b = read_file(path);
    struct solution_file *solution = parse_solution_byte_string(b);
    if (!solution)
        free(b.bytes);
    solution->owns_bytes = true;
    return solution;
}

void free_solution_file(struct solution_file *solution)
{
    for (uint32_t i = 0; i < solution->number_of_parts; ++i) {
        free(solution->parts[i].instructions);
        free(solution->parts[i].track_hexes);
        free(solution->parts[i].conduit_hexes);
    }
    free(solution->parts);
    if (solution->owns_bytes)
        free(solution->bytes);
    free(solution);
}
