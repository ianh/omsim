#include "area.h"

#include "sim.h"
#include <math.h>

// xx remove
#include <assert.h>
#include <stdio.h>

static const double ATOM_RADIUS = 29./82;
static const double SQRT3_2 = 0.8660254037844386;

// xx share this with collision code when that exists
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

static bool intersection(double swing_radius_squared, struct xy_vector center, struct xy_vector *hit)
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

void record_swing_area(struct board *board, struct vector position, struct vector base, int rotation)
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
    // printf("swing radius: %f (%d %d)\n", r, p.u, p.v);
    while (true) {
        // convert back to grid coordinates.
        int32_t cell_u = (int32_t)floor(current.y / SQRT3_2);
        int32_t cell_v = (int32_t)floor(current.x - 0.5 * current.y / SQRT3_2);
        struct xy_vector min = { 0 };
        struct vector min_cell = { 0 };
        // find the intersection point for the hex at each corner of the
        // current grid cell, using the ccw() function to determine which
        // hex the swing arc will encounter next.
        for (int32_t u = cell_u - 1; u <= cell_u + 2; ++u) {
            for (int32_t v = cell_v - 1; v <= cell_v + 2; ++v) {
                struct xy_vector center = to_xy((struct vector){ u, v });
                struct xy_vector hit;
                if (!intersection(r2, center, &hit))
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

void replay_swing_area(struct board *board)
{
    // struct swing_area_table *table = &board->swing_area_table;
    // for (uint32_t i = 0; i < table->capacity; ++i) {
    //     struct swing_area_subtable *subtable = &table->subtables[i];
    //     if (vectors_equal(subtable->offset, zero_vector))
    //         continue;
    //     struct vector p = subtable->offset;
    //     struct xy_vector current = to_xy(p);
    //     struct xy_vector end = to_xy((struct vector){ p.u + p.v, -p.u });
    //     double r2 = xy_len2(current);
    //     double r = sqrt(r2);
    //     while (ccw(current, end) > 0) {
    //         // convert back to grid coordinates.
    //         int32_t cell_u = (int32_t)floor(current.y / SQRT3_2);
    //         int32_t cell_v = (int32_t)floor(current.x - 0.5 * current.y / SQRT3_2);
    //         struct xy_vector min = { 0 };
    //         struct vector min_cell = { 0 };
    //         // find the intersection point for the hex at each corner of the
    //         // current grid cell, using the ccw() function to determine which
    //         // hex the swing arc will encounter next.
    //         for (int32_t u = cell_u; u <= cell_u + 1; ++u) {
    //             for (int32_t v = cell_v; v <= cell_v + 1; ++v) {
    //                 struct xy_vector hit;
    //                 if (!intersection(r, u, v, &hit))
    //                     continue;
    //                 if (ccw(current, hit) <= 0)
    //                     continue;
    //                 if ((min.x != 0 || min.y != 0) && ccw(hit, min) <= 0)
    //                     continue;
    //                 min = hit;
    //                 min_cell.u = u;
    //                 min_cell.v = v;
    //             }
    //         }
    //         for (uint32_t j = 0; j < subtable->capacity; ++j) {
    //             if (vectors_equal(subtable->swings[j].direction, zero_vector))
    //                 continue;
    //             // xx compress the table
    //             struct vector base = subtable->swings[j].base;
    //             struct vector du = subtable->swings[j].direction;
    //             struct vector dv = (struct vector){ du.u + du.v, -du.u };
    //             mark_used_area(board, (struct vector){
    //                 base.u + min_cell.u * du.u + min_cell.v * dv.u,
    //                 base.v + min_cell.u * du.v + min_cell.v * dv.v,
    //             });
    //         }
    //         // advance the current point to the point of encounter.
    //         current = min;
    //     }
    // }
}
