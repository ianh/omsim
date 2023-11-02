#include "collision.h"

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct xy_vector {
    float x;
    float y;
};

struct collider {
    struct xy_vector center;
    float radius;
};

struct collider_list {
    struct collider *colliders;
    size_t cursor;
    size_t length;

    bool collision;
    struct vector *collision_location;
};

// based on image dimensions in textures/board/hex_tile.lighting/
static const float hexSizeX = 82;
static const float hexSizeY = 71;

static const float atomRadius = 29;
static const float producedAtomRadius = 15;
static const float armBaseRadius = 20;

static struct xy_vector to_xy(struct vector p)
{
    return (struct xy_vector){
        hexSizeX * (p.u + 0.5f * p.v),
        hexSizeY * p.v,
    };
}
static struct vector from_xy(struct xy_vector xy)
{
    float u = xy.x / hexSizeX - 0.5f * xy.y / hexSizeY;
    float v = xy.y / hexSizeY;
    float w = 0.f - u - v;
    int ui = (int)round(u);
    int vi = (int)round(v);
    int wi = (int)round(w);
    // the original code uses double-precision fabs, but since the argument is
    // single-precision, fabsf should return the same result.  it's faster not
    // to have to convert to double and back.
    float uf = fabsf(u - (float)ui);
    float vf = fabsf(v - (float)vi);
    float wf = fabsf(w - (float)wi);
    if (uf > vf && uf > wf)
        ui = -vi - wi;
    else if (vf > wf)
        vi = -ui - wi;
    return (struct vector){ ui, vi };
}
static float xy_len2(struct xy_vector xy)
{
    return xy.x * xy.x + xy.y * xy.y;
}
static struct xy_vector xy_sub(struct xy_vector a, struct xy_vector b)
{
    return (struct xy_vector){ a.x - b.x, a.y - b.y };
}
static float xy_len(struct xy_vector xy)
{
    return (float)sqrt(xy_len2(xy));
}
static float xy_dist(struct xy_vector a, struct xy_vector b)
{
    return xy_len(xy_sub(a, b));
}
static float to_radians(int32_t r)
{
    return (float)(r * 60) * ((float)M_PI / 180.f);
}

static void mark_area_and_check_board(struct collider_list *list, struct board *board, struct collider collider, int32_t u, int32_t v)
{
    struct vector p = { u, v };
    struct xy_vector center = to_xy(p);
    float dist = xy_dist(collider.center, center);
    if (!(dist < collider.radius + atomRadius))
        return;
    // mark area.  also, this *could* be a collision.
    atom a = mark_used_area(board, p);
    if (a & REMOVED)
        return;
    if (a & BEING_PRODUCED) {
        // atoms have a somewhat smaller collision radius as they emerge from a glyph.
        if (!(dist < collider.radius + producedAtomRadius))
            return;
    }
    list->collision = true;
    if (list->collision_location)
        *list->collision_location = p;
}

__attribute__((always_inline))
static void add_collider_inline(struct collider_list *list, struct board *board, struct collider collider)
{
    struct vector p = from_xy(collider.center);
    mark_area_and_check_board(list, board, collider, p.u, p.v);
    mark_area_and_check_board(list, board, collider, p.u + 1, p.v);
    mark_area_and_check_board(list, board, collider, p.u, p.v + 1);
    mark_area_and_check_board(list, board, collider, p.u - 1, p.v + 1);
    mark_area_and_check_board(list, board, collider, p.u - 1, p.v);
    mark_area_and_check_board(list, board, collider, p.u, p.v - 1);
    mark_area_and_check_board(list, board, collider, p.u + 1, p.v - 1);
    list->colliders[list->length++] = collider;
    for (size_t i = 0; i < list->cursor; ++i) {
        struct collider other = list->colliders[i];
        if (!(xy_dist(other.center, collider.center) < other.radius + collider.radius))
            continue;
        list->collision = true;
        if (list->collision_location)
            *list->collision_location = p;
        break;
    }
}

static void add_collider(struct collider_list *list, struct board *board, struct collider collider)
{
    add_collider_inline(list, board, collider);
}

bool collision(struct solution *solution, struct board *board, float increment, struct vector *collision_location)
{
    size_t number_of_colliders = solution->number_of_arms;
    if (!board->ignore_swing_area)
        number_of_colliders += board->moving_atoms.length;
    struct collider_list list = {
        .colliders = calloc(number_of_colliders, sizeof(struct collider)),
        .collision_location = collision_location,
    };
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        if (!vectors_equal(solution->arms[i].movement, zero_vector))
            continue;
        add_collider(&list, board, (struct collider){
            .center = to_xy(solution->arms[i].position),
            .radius = armBaseRadius,
        });
        list.cursor = list.length;
    }
    list.cursor = list.length;
    size_t fixed_colliders = list.length;
    for (float progress = increment; progress < 1.f; progress += increment) {
        list.cursor = fixed_colliders;
        list.length = fixed_colliders;
        for (size_t i = 0; i < solution->number_of_arms; ++i) {
            if (vectors_equal(solution->arms[i].movement, zero_vector))
                continue;
            struct xy_vector xy = to_xy(solution->arms[i].position);
            struct xy_vector tr = to_xy(solution->arms[i].movement);
            xy.x -= tr.x * (1.f - progress);
            xy.y -= tr.y * (1.f - progress);
            add_collider(&list, board, (struct collider){
                .center = xy,
                .radius = armBaseRadius,
            });
            list.cursor = list.length;
        }
        // xx rate should still check swing collision/area, but in a different
        // way that isn't so slow...
        if (board->ignore_swing_area)
            continue;
        size_t atom_index = 0;
        for (size_t i = 0; i < board->movements.length; ++i) {
            struct movement m = board->movements.movements[i];
            struct xy_vector v = to_xy(m.base);
            float rotation = to_radians(m.rotation);
            float armRotation = 0;
            float r = (0.f - rotation) * (1.f - progress);
            switch (m.type & 3) {
            case SWING_MOVEMENT:
                armRotation = rotation;
                break;
            case TRACK_MOVEMENT: {
                struct xy_vector tr = to_xy(m.translation);
                v.x -= tr.x * (1.f - progress);
                v.y -= tr.y * (1.f - progress);
                break;
            }
            default:
                break;
            }
            float baseRotation = to_radians(m.base_rotation);
            struct xy_vector g = to_xy(m.grabber_offset);
            if (m.type & IS_PISTON) {
                float len = xy_len(g);
                float newlen = xy_len(g) - ((float)m.piston_extension) * hexSizeX * (1.f - progress);
                g.x = newlen / len * g.x;
                g.y = newlen / len * g.y;
            }
            float grabberRotation = baseRotation - armRotation * (1.f - progress);
            float grx = (float)cos(grabberRotation);
            float gry = (float)sin(grabberRotation);
            v.x += g.x * grx - g.y * gry;
            v.y += g.x * gry + g.y * grx;
            float rx = (float)cos(r);
            float ry = (float)sin(r);
            for (size_t j = 0; j < m.number_of_atoms; ++j) {
                struct vector p = board->moving_atoms.atoms_at_positions[atom_index++].position;
                p.u -= m.absolute_grab_position.u;
                p.v -= m.absolute_grab_position.v;
                struct xy_vector xy = to_xy(p);
                add_collider_inline(&list, board, (struct collider){
                    .center = {
                        v.x + (xy.x * rx - xy.y * ry),
                        v.y + (xy.x * ry + xy.y * rx),
                    },
                    .radius = atomRadius,
                });
            }
            list.cursor = list.length;
        }

        // printf("[");
        // for (size_t i = 0; i < list.length; ++i) {
        //     if (i > 0)
        //         printf(",");
        //     printf("[%f,%f,%f]", list.colliders[i].radius, list.colliders[i].center.x, list.colliders[i].center.y);
        // }
        // printf("],");
    }

    free(list.colliders);
    return list.collision;
}
