#include "decode.h"
#include "parse.h"
#include "sim.h"
#include <limits.h>
#include <stdlib.h>

static unsigned char stabilized_water[] = {
  0x03, 0x00, 0x00, 0x00, 0x10, 0x53, 0x54, 0x41, 0x42, 0x49, 0x4c, 0x49,
  0x5a, 0x45, 0x44, 0x20, 0x57, 0x41, 0x54, 0x45, 0x52, 0xce, 0xf9, 0x01,
  0x02, 0x01, 0x00, 0x10, 0x01, 0x0f, 0x17, 0xc0, 0x07, 0x00, 0x00, 0x00,
  0x00, 0x02, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x05, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x01,
  0x00, 0x00, 0x05, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00,
  0x01, 0x00, 0x01, 0x00, 0x00, 0x00, 0x00
};

static unsigned char armor_filament[] = {
  0x03, 0x00, 0x00, 0x00, 0x0e, 0x41, 0x52, 0x4d, 0x4f, 0x52, 0x20, 0x46,
  0x49, 0x4c, 0x41, 0x4d, 0x45, 0x4e, 0x54, 0xce, 0xf9, 0x01, 0x02, 0x01,
  0x00, 0x10, 0x01, 0x0f, 0x57, 0xc0, 0x07, 0x00, 0x00, 0x00, 0x00, 0x02,
  0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x06, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x01, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x0a, 0x00, 0x00,
  0x0a, 0x01, 0x00, 0x0f, 0x02, 0x00, 0x02, 0x00, 0x00, 0x00, 0x01, 0x01,
  0x00, 0x02, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x01, 0x00, 0x00, 0x00,
  0x00
};

static struct solution_file *parse_alt_solution_byte_string(struct byte_string b);

static struct puzzle_file *pf;

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t length)
{
    if (!pf) {
        // struct byte_string b = { stabilized_water, sizeof(stabilized_water) };
        struct byte_string b = { armor_filament, sizeof(armor_filament) };
        pf = parse_puzzle_byte_string(b);
        if (!pf)
            abort();
    }
    struct byte_string input = { (unsigned char *)data, length };
    struct solution_file *sf = parse_alt_solution_byte_string(input);
    if (!sf)
        return 0;
    struct solution solution = { 0 };
    struct board board = { 0 };
    if (decode_solution(&solution, pf, sf, 0)) {
        initial_setup(&solution, &board, sf->area);
        while (board.cycle < 200 && !board.complete) {
            // printf("-- %llu %u %u\n", board.cycle, board.capacity, board.used);
            cycle(&solution, &board);
            if (board.collision)
                break;
        }
    }
    destroy(&solution, &board);
    free_solution_file(sf);
    return 0;
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

static struct byte_string part_name_for_byte(int byte)
{
    switch (byte) {
    case 0:
        return (struct byte_string){ (unsigned char *)"arm1", strlen("arm1") };
    case 1:
        return (struct byte_string){ (unsigned char *)"arm2", strlen("arm2") };
    case 2:
        return (struct byte_string){ (unsigned char *)"arm3", strlen("arm3") };
    case 3:
        return (struct byte_string){ (unsigned char *)"arm6", strlen("arm6") };
    case 4:
        return (struct byte_string){ (unsigned char *)"piston", strlen("piston") };
    case 5:
        return (struct byte_string){ (unsigned char *)"baron", strlen("baron") };
    case 6:
        return (struct byte_string){ (unsigned char *)"track", strlen("track") };
    case 7:
        return (struct byte_string){ (unsigned char *)"bonder", strlen("bonder") };
    case 8:
        return (struct byte_string){ (unsigned char *)"unbonder", strlen("unbonder") };
    case 9:
        return (struct byte_string){ (unsigned char *)"bonder-prisma", strlen("bonder-prisma") };
    case 10:
        return (struct byte_string){ (unsigned char *)"bonder-speed", strlen("bonder-speed") };
    case 11:
        return (struct byte_string){ (unsigned char *)"glyph-calcification", strlen("glyph-calcification") };
    case 12:
        return (struct byte_string){ (unsigned char *)"glyph-dispersion", strlen("glyph-dispersion") };
    case 13:
        return (struct byte_string){ (unsigned char *)"glyph-disposal", strlen("glyph-disposal") };
    case 14:
        return (struct byte_string){ (unsigned char *)"glyph-duplication", strlen("glyph-duplication") };
    case 15:
        return (struct byte_string){ (unsigned char *)"glyph-life-and-death", strlen("glyph-life-and-death") };
    case 16:
        return (struct byte_string){ (unsigned char *)"glyph-marker", strlen("glyph-marker") };
    case 17:
        return (struct byte_string){ (unsigned char *)"glyph-projection", strlen("glyph-projection") };
    case 18:
        return (struct byte_string){ (unsigned char *)"glyph-purification", strlen("glyph-purification") };
    case 19:
        return (struct byte_string){ (unsigned char *)"glyph-unification", strlen("glyph-unification") };
    case 20:
        return (struct byte_string){ (unsigned char *)"input", strlen("input") };
    case 21:
        return (struct byte_string){ (unsigned char *)"out-rep", strlen("out-rep") };
    case 22:
        return (struct byte_string){ (unsigned char *)"out-std", strlen("out-std") };
    case 23:
        return (struct byte_string){ (unsigned char *)"pipe", strlen("pipe") };
    default:
        return (struct byte_string){ 0 };
    }
}

static struct solution_file *parse_alt_solution_byte_string(struct byte_string b)
{
    struct solution_file *solution = calloc(1, sizeof(struct solution_file));
    solution->bytes = b.bytes;
    solution->parts = calloc(999, sizeof(struct solution_part));
    while (b.length > 0 && solution->number_of_parts < 999) {
        uint32_t n = solution->number_of_parts++;
        struct solution_part *part = &solution->parts[n];
        int part_name_byte = read_byte(&b);
        part->name = part_name_for_byte(part_name_byte);
        part->position[0] = read_byte(&b);
        part->position[1] = read_byte(&b);
        part->rotation = read_byte(&b);
        if (part_name_byte <= 5) {
            part->size = read_byte(&b);
            if (part->size > 3)
                part->size = 3;
            uint32_t instruction_offset = read_byte(&b);
            part->instructions = calloc(999, sizeof(struct solution_instruction));
            while (part->number_of_instructions < 999) {
                uint32_t m = part->number_of_instructions;
                struct solution_instruction *inst = &part->instructions[m];
                inst->index = instruction_offset + m;
                inst->instruction = read_byte(&b);
                if (inst->instruction == 0)
                    break;
                part->number_of_instructions++;
            }
        } else if (part_name_byte == 6) {
            part->track_hexes = calloc(99, sizeof(struct solution_hex_offset));
            while (part->number_of_track_hexes < 99) {
                uint32_t m = part->number_of_track_hexes++;
                struct solution_hex_offset *hex = &part->track_hexes[m];
                hex->offset[0] = read_byte(&b);
                hex->offset[1] = read_byte(&b);
                if (read_byte(&b) < 100)
                    break;
            }
        } else if (part_name_byte >= 20 && part_name_byte <= 22)
            part->which_input_or_output = read_byte(&b);
        else if (part_name_byte == 23) {
            part->conduit_id = read_byte(&b);
            part->conduit_hexes = calloc(99, sizeof(struct solution_hex_offset));
            while (part->number_of_conduit_hexes < 99) {
                uint32_t m = part->number_of_conduit_hexes++;
                struct solution_hex_offset *hex = &part->conduit_hexes[m];
                hex->offset[0] = read_byte(&b);
                hex->offset[1] = read_byte(&b);
                if (read_byte(&b) < 100)
                    break;
            }
        }
    }
    return solution;
}
