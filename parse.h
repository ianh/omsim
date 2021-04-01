#ifndef OM_PARSE_H
#define OM_PARSE_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

struct byte_string {
    unsigned char *bytes;
    size_t length;
};

struct puzzle_file;
struct solution_file;

struct puzzle_file *parse_puzzle_file(const char *path);
struct puzzle_file *parse_puzzle_byte_string(struct byte_string b);
void free_puzzle_file(struct puzzle_file *puzzle);

struct solution_file *parse_solution_file(const char *path);
struct solution_file *parse_solution_byte_string(struct byte_string b);
void free_solution_file(struct solution_file *solution);

static inline bool byte_string_is(struct byte_string b, const char *str)
{
    size_t len = strlen(str);
    return b.length == len && !memcmp(b.bytes, str, len);
}

struct puzzle_atom;
struct puzzle_bond;
struct puzzle_molecule;
struct puzzle_file {
    struct byte_string name;

    uint64_t creator;
    uint64_t parts_available;

    uint32_t number_of_inputs;
    struct puzzle_molecule *inputs;

    uint32_t number_of_outputs;
    struct puzzle_molecule *outputs;

    // xx name of this field?
    uint32_t output_scale;

    // xx cabinets

    void *bytes;
    bool owns_bytes;
};
struct puzzle_molecule {
    uint32_t number_of_atoms;
    struct puzzle_atom *atoms;

    uint32_t number_of_bonds;
    struct puzzle_bond *bonds;
};
struct puzzle_atom {
    unsigned char type;
    signed char offset[2];
};
struct puzzle_bond {
    unsigned char type;
    signed char from[2];
    signed char to[2];
};

struct solution_hex_offset;
struct solution_instruction;
struct solution_part;
struct solution_file {
    struct byte_string puzzle;
    struct byte_string name;

    uint32_t cycles;
    uint32_t cost;
    uint32_t area;
    uint32_t instructions;

    uint32_t number_of_parts;
    struct solution_part *parts;

    void *bytes;
    bool owns_bytes;
};
struct solution_part {
    struct byte_string name;

    int32_t position[2];
    int32_t size;
    int32_t rotation;

    // used for part names "input", "out-std", and "out-rep".
    uint32_t which_input_or_output;

    uint32_t number_of_instructions;
    struct solution_instruction *instructions;

    // nonzero only when the part name is "track".
    uint32_t number_of_track_hexes;
    struct solution_hex_offset *track_hexes;

    uint32_t arm_number;

    // nonzero only when the part name is "pipe".
    uint32_t conduit_id;
    uint32_t number_of_conduit_hexes;
    struct solution_hex_offset *conduit_hexes;
};
struct solution_instruction {
    int32_t index;
    char instruction;
};
struct solution_hex_offset {
    int32_t offset[2];
};


#endif
