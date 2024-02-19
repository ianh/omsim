#include "collision.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

// M_PI isn't standard, so we define it ourselves here.
#ifdef M_PI
#undef M_PI
#endif
#define M_PI 3.14159265358979323846

// only try to find collisions up to this number of periods.
#define CHAIN_ATOM_PERIOD_LIMIT 500

struct xy_vector {
    float x;
    float y;
};

struct xy_rect {
    float min_x;
    float max_x;
    float min_y;
    float max_y;
};
static const struct xy_rect empty_rect = { INFINITY, -INFINITY, INFINITY, -INFINITY };

struct collider {
    struct xy_vector center;
    float radius;
};
struct chain_atom_collider {
    struct xy_vector v;
    struct xy_vector r;
    struct vector offset_from_grab;
    struct vector per_period_motion;
    size_t movement_index;
    int32_t extra_periods;
    uint32_t area_direction;
    bool in_repeating_segment;
};

struct collider_list {
    struct collider *colliders;
    size_t cursor;
    size_t length;

    struct chain_atom_collider *chain_atom_colliders;
    size_t number_of_chain_atom_colliders;

    struct xy_rect bounding_box;
    struct xy_rect bounding_box_up_to_cursor;

    bool collision;
    struct vector *collision_location;

    bool ignore_board;
    bool bounding_box_unused;
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

static struct xy_rect xy_rect_add_collider(struct xy_rect rect, struct collider collider)
{
    if (collider.center.x - collider.radius < rect.min_x)
        rect.min_x = collider.center.x - collider.radius;
    if (collider.center.x + collider.radius > rect.max_x)
        rect.max_x = collider.center.x + collider.radius;
    if (collider.center.y - collider.radius < rect.min_y)
        rect.min_y = collider.center.y - collider.radius;
    if (collider.center.y + collider.radius > rect.max_y)
        rect.max_y = collider.center.y + collider.radius;
    return rect;
}
static bool xy_rect_contains_collider(struct xy_rect rect, struct collider collider)
{
    if (collider.center.x + collider.radius < rect.min_x)
        return false;
    if (collider.center.x - collider.radius > rect.max_x)
        return false;
    if (collider.center.y + collider.radius < rect.min_y)
        return false;
    if (collider.center.y - collider.radius > rect.max_y)
        return false;
    return true;
}

static bool xy_rect_intersects_ray(struct xy_rect rect, struct xy_vector origin, double motion_delta_x, double motion_delta_y, float radius)
{
    if (motion_delta_x < 0) {
        float tmp = rect.min_x;
        rect.min_x = -rect.max_x;
        rect.max_x = -tmp;
        motion_delta_x = -motion_delta_x;
    }
    if (motion_delta_y < 0) {
        float tmp = rect.min_y;
        rect.min_y = -rect.max_y;
        rect.max_y = -tmp;
        motion_delta_y = -motion_delta_y;
    }
    double max_x = rect.max_x + radius - origin.x;
    double max_y = rect.max_y + radius - origin.y;
    if (max_x < 0 || max_y < 0)
        return false;
    double min_x = rect.min_x - radius - origin.x;
    double min_y = rect.min_y - radius - origin.y;
    if (min_x * motion_delta_y > max_y * motion_delta_x || min_y * motion_delta_x > max_x * motion_delta_y)
        return false;
    return true;
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
static void add_collider(struct collider_list *list, struct board *board, struct collider collider)
{
    if (!list->ignore_board) {
        struct vector p = from_xy(collider.center);
        mark_area_and_check_board(list, board, collider, p.u, p.v);
        mark_area_and_check_board(list, board, collider, p.u + 1, p.v);
        mark_area_and_check_board(list, board, collider, p.u, p.v + 1);
        mark_area_and_check_board(list, board, collider, p.u - 1, p.v + 1);
        mark_area_and_check_board(list, board, collider, p.u - 1, p.v);
        mark_area_and_check_board(list, board, collider, p.u, p.v - 1);
        mark_area_and_check_board(list, board, collider, p.u + 1, p.v - 1);
    }
    list->colliders[list->length++] = collider;
    if (!list->bounding_box_unused)
        list->bounding_box = xy_rect_add_collider(list->bounding_box, collider);
    if (!xy_rect_contains_collider(list->bounding_box_up_to_cursor, collider))
        return;
    for (size_t i = 0; i < list->cursor; ++i) {
        struct collider other = list->colliders[i];
        if (!(xy_dist(other.center, collider.center) < other.radius + collider.radius))
            continue;
        list->collision = true;
        struct vector p = from_xy(collider.center);
        if (list->collision_location)
            *list->collision_location = p;
        break;
    }
}

static void add_chain_atom_collider(struct collider_list *list, struct board *board, size_t movement_index, uint32_t chain_atom)
{
    struct chain_atom ca = board->chain_atoms[chain_atom];
    struct vector absolute_grab_position = zero_vector;
    if (movement_index < board->movements.length)
        absolute_grab_position = board->movements.movements[movement_index].absolute_grab_position;
    list->chain_atom_colliders[list->number_of_chain_atom_colliders++] = (struct chain_atom_collider){
        .r = { 1, 0 },
        .offset_from_grab = { ca.current_position.u - absolute_grab_position.u, ca.current_position.v - absolute_grab_position.v },
        .per_period_motion = { ca.current_position.u - ca.original_position.u, ca.current_position.v - ca.original_position.v },
        .movement_index = movement_index,
        .area_direction = ca.area_direction,
        .in_repeating_segment = ca.flags & CHAIN_ATOM_IN_REPEATING_SEGMENT,
    };
}

static struct xy_vector chain_atom_center_for_period(struct chain_atom_collider collider, int32_t period)
{
    struct xy_vector xy = to_xy((struct vector){
        collider.offset_from_grab.u + period * collider.per_period_motion.u,
        collider.offset_from_grab.v + period * collider.per_period_motion.v,
    });
    return (struct xy_vector){
        collider.v.x + (xy.x * collider.r.x - xy.y * collider.r.y),
        collider.v.y + (xy.x * collider.r.y + xy.y * collider.r.x),
    };
}

static int compare_vectors(struct vector a, struct vector b)
{
    if (a.u < b.u)
        return -1;
    if (a.u > b.u)
        return 1;
    if (a.v < b.v)
        return -1;
    if (a.v > b.v)
        return 1;
    return 0;
}

static void get_period_and_equivalence_class(struct chain_atom_collider collider, int32_t *period, struct vector *equivalence_class)
{
    if (collider.per_period_motion.u != 0)
        *period = collider.offset_from_grab.u / collider.per_period_motion.u;
    else if (collider.per_period_motion.v != 0)
        *period = collider.offset_from_grab.v / collider.per_period_motion.v;
    else
        *period = 0;
    equivalence_class->u = collider.offset_from_grab.u - *period * collider.per_period_motion.u;
    equivalence_class->v = collider.offset_from_grab.v - *period * collider.per_period_motion.v;
}

static int compare_chain_atom_colliders_modulo_period(struct chain_atom_collider a, struct chain_atom_collider b, int32_t *a_period, int32_t *b_period)
{
    int cmp = compare_vectors(a.per_period_motion, b.per_period_motion);
    if (cmp != 0)
        return cmp;
    struct vector a_equivalence_class, b_equivalence_class;
    get_period_and_equivalence_class(a, a_period, &a_equivalence_class);
    get_period_and_equivalence_class(b, b_period, &b_equivalence_class);
    return compare_vectors(a_equivalence_class, b_equivalence_class);
}

static int compare_chain_atom_colliders_within_movement(const void *a, const void *b)
{
    int32_t a_period, b_period;
    int cmp = compare_chain_atom_colliders_modulo_period(*(struct chain_atom_collider *)a, *(struct chain_atom_collider *)b, &a_period, &b_period);
    if (cmp != 0)
        return cmp;
    return a_period - b_period;
}

static void combine_chain_atoms(struct collider_list *list, size_t start_index)
{
    qsort(list->chain_atom_colliders + start_index, list->number_of_chain_atom_colliders - start_index, sizeof(struct chain_atom_collider), compare_chain_atom_colliders_within_movement);
    size_t removed = 0;
    size_t combine_into = start_index;
    for (size_t i = start_index + 1; i < list->number_of_chain_atom_colliders; ++i) {
        int32_t extra_periods = list->chain_atom_colliders[combine_into].extra_periods;
        int32_t prev_period, period;
        int cmp = compare_chain_atom_colliders_modulo_period(list->chain_atom_colliders[combine_into], list->chain_atom_colliders[i], &prev_period, &period);
        if (cmp == 0 && period - prev_period == extra_periods + 1 && list->chain_atom_colliders[combine_into].in_repeating_segment) {
            list->chain_atom_colliders[combine_into].extra_periods++;
            removed++;
        } else {
            list->chain_atom_colliders[i - removed] = list->chain_atom_colliders[i];
            combine_into = i - removed;
        }
    }
    list->number_of_chain_atom_colliders -= removed;
}

static void get_per_period_motion(struct chain_atom_collider collider, double *per_period_motion_x, double *per_period_motion_y)
{
    double x = (double)hexSizeX * (collider.per_period_motion.u + 0.5 * collider.per_period_motion.v);
    double y = (double)hexSizeY * collider.per_period_motion.v;
    *per_period_motion_x = x * collider.r.x - y * collider.r.y;
    *per_period_motion_y = x * collider.r.y + y * collider.r.x;
}

static int32_t minimum_approach_period(double motion_delta_x, double motion_delta_y, struct xy_vector origin_a, struct xy_vector origin_b)
{
    double period = motion_delta_x * (origin_b.x - origin_a.x) + motion_delta_y * (origin_b.y - origin_a.y);
    period /= motion_delta_x * motion_delta_x + motion_delta_y * motion_delta_y;
    if (!(period >= 0))
        period = 0;
    if (period > CHAIN_ATOM_PERIOD_LIMIT)
        period = CHAIN_ATOM_PERIOD_LIMIT;
    return lrint(period);
}

// static void print_collider(struct chain_atom_collider a)
// {
//     printf("v=(%f %f) r=(%f %f) ofg=(%d %d) ppm=(%d %d) ep=%d irs=%d\n", a.v.x, a.v.y, a.r.x, a.r.y,
//         a.offset_from_grab.u, a.offset_from_grab.v,
//         a.per_period_motion.u, a.per_period_motion.v,
//         a.extra_periods, a.in_repeating_segment);
// }

__attribute__((noinline))
static void resolve_chain_atom_collisions(struct collider_list *list)
{
    for (size_t i = 0; i < list->number_of_chain_atom_colliders; ++i) {
        struct chain_atom_collider a = list->chain_atom_colliders[i];
        struct xy_vector origin_a = chain_atom_center_for_period(a, 0);
        double motion_ax = 0;
        double motion_ay = 0;
        get_per_period_motion(a, &motion_ax, &motion_ay);
        double motion_length_a = sqrt(motion_ax * motion_ax + motion_ay * motion_ay);
        double perp_ax = -motion_ay / motion_length_a * atomRadius * 2;
        double perp_ay = motion_ax / motion_length_a * atomRadius * 2;
        // first, check for collisions between the chain atom and normal colliders.
        if (xy_rect_intersects_ray(list->bounding_box, origin_a, motion_ax, motion_ay, atomRadius)) {
            for (size_t j = 0; j < list->length; ++j) {
                struct collider b = list->colliders[j];
                int32_t period = minimum_approach_period(motion_ax, motion_ay, origin_a, b.center);
                if (xy_dist(b.center, chain_atom_center_for_period(a, period)) < b.radius + atomRadius) {
                    list->collision = true;
                    if (list->collision_location)
                        *list->collision_location = from_xy(b.center);
                }
            }
        }
        // then, check for collisions between the chain atom and other chain atoms.
        for (size_t j = 0; j < i; ++j) {
            struct chain_atom_collider b = list->chain_atom_colliders[j];
            if (b.movement_index == a.movement_index)
                break;
            struct xy_vector origin_b = chain_atom_center_for_period(b, 0);
            double motion_bx = 0;
            double motion_by = 0;
            get_per_period_motion(b, &motion_bx, &motion_by);
            if (!a.in_repeating_segment && !b.in_repeating_segment) {
                // if neither atom is in a repeating segment, subtracting their motions reduces this to the stationary case.
                int32_t period = minimum_approach_period(motion_ax - motion_bx, motion_ay - motion_by, origin_a, origin_b);
                if (xy_dist(chain_atom_center_for_period(a, period), chain_atom_center_for_period(b, period)) < atomRadius + atomRadius) {
                    list->collision = true;
                    if (list->collision_location)
                        *list->collision_location = from_xy(chain_atom_center_for_period(a, period));
                }
            } else {
                // one of the two atoms is in a repeating segment.  find where the two lines of motion get close enough to interact.
                int32_t start = 0;
                int32_t end = CHAIN_ATOM_PERIOD_LIMIT;
                double d = motion_bx * motion_ay - motion_ax * motion_by;
                if (fabs(d) < 1e-6) {
                    // the lines of motion are approximately parallel.
                    // skip checking if they're further away than 2 atom radiuses.
                    double closest_period = motion_ax * (origin_b.x - origin_a.x) + motion_ay * (origin_b.y - origin_a.y);
                    closest_period /= motion_ax * motion_ax + motion_ay * motion_ay;
                    double closest_x = origin_b.x - origin_a.x - closest_period * motion_ax;
                    double closest_y = origin_b.y - origin_a.y - closest_period * motion_ay;
                    if (closest_x * closest_x + closest_y * closest_y > 4 * atomRadius * atomRadius)
                        continue;
                } else {
                    // determine the range of periods where intersections can occur.
                    double t1 = ((origin_b.y - origin_a.y - perp_ay) * motion_ax - (origin_b.x - origin_a.x - perp_ax) * motion_ay) / d;
                    double t2 = ((origin_b.y - origin_a.y + perp_ay) * motion_ax - (origin_b.x - origin_a.x + perp_ax) * motion_ay) / d;
                    if (t2 < t1) {
                        double tmp = t1;
                        t1 = t2;
                        t2 = tmp;
                    }
                    if (t1 > CHAIN_ATOM_PERIOD_LIMIT || t2 < 0)
                        continue;
                    if (!(t1 >= 0))
                        t1 = 0;
                    if (!(t2 <= CHAIN_ATOM_PERIOD_LIMIT))
                        t2 = CHAIN_ATOM_PERIOD_LIMIT;
                    start = t1;
                    end = t2;
                }
                // search the range of possible intersecting periods for a collision.
                for (int32_t k = start; k <= end; ++k) {
                    struct xy_vector origin_b_at_period = chain_atom_center_for_period(b, k);
                    int32_t period = minimum_approach_period(motion_ax, motion_ay, origin_a, origin_b_at_period);
                    if (a.in_repeating_segment && !b.in_repeating_segment && period > k + a.extra_periods)
                        period = k + a.extra_periods;
                    if (!a.in_repeating_segment && b.in_repeating_segment && period < k - b.extra_periods)
                        period = k - b.extra_periods;
                    if (xy_dist(origin_b_at_period, chain_atom_center_for_period(a, period)) < atomRadius + atomRadius) {
                        list->collision = true;
                        if (list->collision_location)
                            *list->collision_location = from_xy(origin_b_at_period);
                    }
                }
            }
        }
    }
}

// stolen from python
static int32_t floor_div(int32_t a, int32_t b)
{
    if (a == 0)
        return 0;
    else if (a < 0 == b < 0)
        return labs(a) / labs(b);
    else
        return -1 - (labs(a) - 1) / labs(b);
}

static void mark_area_at_infinity_in_direction(struct linear_area_direction *d, struct xy_vector center, int32_t divisions, int32_t u, int32_t v)
{
    struct vector p = { u, v };
    struct xy_vector grid_center = to_xy(p);
    float dist = xy_dist(center, grid_center);
    if (!(dist < atomRadius + atomRadius))
        return;
    int32_t period = 0;
    if (d->direction.u != 0)
        period = floor_div(p.u, d->direction.u);
    else
        period = floor_div(p.v, d->direction.v);
    p.u -= period * d->direction.u;
    p.v -= period * d->direction.v;
    struct atom_at_position *ap = lookup_atom_at_position(&d->footprint_at_infinity, p);
    if (!(ap->atom & VALID) || ((ap->atom >> 1) > divisions)) {
        ap->position = p;
        ap->atom = VALID | ((atom)divisions << 1);
    }
}

static void mark_area_at_infinity(struct board *board, struct chain_atom_collider collider)
{
    struct linear_area_direction *d = &board->area_directions[collider.area_direction];
    int32_t divisions = 1;
    if (collider.per_period_motion.u != 0)
        divisions = d->direction.u / collider.per_period_motion.u;
    else
        divisions = d->direction.v / collider.per_period_motion.v;
    for (int32_t i = 0; i < divisions; ++i) {
        struct xy_vector center = chain_atom_center_for_period(collider, i);
        struct vector p = from_xy(center);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u, p.v);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u + 1, p.v);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u, p.v + 1);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u - 1, p.v + 1);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u - 1, p.v);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u, p.v - 1);
        mark_area_at_infinity_in_direction(d, center, divisions, p.u + 1, p.v - 1);
    }
}

bool collision(struct solution *solution, struct board *board, float increment, struct vector *collision_location)
{
    enum chain_mode chain_mode = board->chain_mode;
    size_t number_of_colliders = solution->number_of_arms + solution->number_of_cabinet_walls + board->moving_atoms.length;
    size_t number_of_chain_atom_colliders = 0;
    if (chain_mode == EXTEND_CHAIN) {
        for (size_t i = 0; i < BOARD_CAPACITY(board); ++i) {
            struct atom_at_position ap = board->grid.atoms_at_positions[i];
            if (!(ap.atom & VALID) || (ap.atom & REMOVED) || (ap.atom & IS_CHAIN_ATOM))
                continue;
            number_of_colliders++;
        }
        for (uint32_t i = 0; i < board->number_of_chain_atoms; ++i) {
            struct chain_atom *ca = &board->chain_atoms[i];
            if (ca->prev_in_list)
                number_of_chain_atom_colliders++;
        }
    }
    struct collider_list list = {
        .colliders = calloc(number_of_colliders, sizeof(struct collider)),
        .chain_atom_colliders = calloc(number_of_chain_atom_colliders, sizeof(struct chain_atom_collider)),
        .bounding_box = empty_rect,
        .bounding_box_up_to_cursor = empty_rect,
        .collision_location = collision_location,
    };
    if (chain_mode == EXTEND_CHAIN) {
        list.ignore_board = true;
        for (size_t i = 0; i < BOARD_CAPACITY(board); ++i) {
            struct atom_at_position ap = board->grid.atoms_at_positions[i];
            if (!(ap.atom & VALID) || (ap.atom & REMOVED))
                continue;
            if (ap.atom & IS_CHAIN_ATOM)
                add_chain_atom_collider(&list, board, SIZE_MAX, lookup_chain_atom(board, ap.position));
            else {
                add_collider(&list, board, (struct collider){
                    .center = to_xy(ap.position),
                    .radius = (ap.atom & BEING_PRODUCED) ? producedAtomRadius : atomRadius,
                });
            }
        }
        list.ignore_board = board->area_growth_order >= GROWTH_LINEAR;
        combine_chain_atoms(&list, 0);
        list.cursor = list.length;
        list.bounding_box_up_to_cursor = list.bounding_box;
    }
    size_t fixed_chain_atom_colliders = list.number_of_chain_atom_colliders;
    size_t atom_index = 0;
    struct xy_vector *moving_atom_offsets = calloc(board->moving_atoms.length, sizeof(struct xy_vector));
    for (size_t i = 0; i < board->movements.length; ++i) {
        struct movement m = board->movements.movements[i];
        for (size_t j = 0; j < m.number_of_atoms; ++j) {
            if (chain_mode == EXTEND_CHAIN && board->moving_atoms.atoms_at_positions[atom_index].atom & IS_CHAIN_ATOM)
                continue;
            struct vector p = board->moving_atoms.atoms_at_positions[atom_index].position;
            p.u -= m.absolute_grab_position.u;
            p.v -= m.absolute_grab_position.v;
            moving_atom_offsets[atom_index++] = to_xy(p);
        }
        if (chain_mode == EXTEND_CHAIN) {
            size_t first_collider_index = list.number_of_chain_atom_colliders;
            uint32_t chain = m.first_chain_atom;
            while (chain != UINT32_MAX) {
                add_chain_atom_collider(&list, board, i, chain);
                chain = board->chain_atoms[chain].next_in_list;
            }
            combine_chain_atoms(&list, first_collider_index);
        }
    }
    for (size_t i = 0; i < solution->number_of_arms; ++i) {
        if (!vectors_equal(solution->arms[i].movement, zero_vector))
            continue;
        add_collider(&list, board, (struct collider){
            .center = to_xy(solution->arms[i].position),
            .radius = armBaseRadius,
        });
        list.cursor = list.length;
        list.bounding_box_up_to_cursor = list.bounding_box;
    }
    for (size_t i = 0; i < solution->number_of_cabinet_walls; ++i) {
        add_collider(&list, board, (struct collider){
            .center = to_xy(solution->cabinet_walls[i]),
            .radius = armBaseRadius,
        });
    }
    size_t fixed_colliders = list.length;
    struct xy_rect fixed_bounding_box = list.bounding_box;
    for (float progress = increment; progress < 1.f; progress += increment) {
        list.bounding_box_unused = false;
        list.cursor = fixed_colliders;
        list.length = fixed_colliders;
        list.bounding_box = fixed_bounding_box;
        list.bounding_box_up_to_cursor = fixed_bounding_box;
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
            list.bounding_box_up_to_cursor = list.bounding_box;
        }
        size_t atom_index = 0;
        size_t chain_atom_cursor = fixed_chain_atom_colliders;
        for (size_t i = 0; i < board->movements.length; ++i) {
            // skip updating the bounding box for the last movement.
            list.bounding_box_unused = i + 1 == board->movements.length;
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
                // xx does this have a perf impact?
                if (chain_mode == EXTEND_CHAIN && board->moving_atoms.atoms_at_positions[atom_index].atom & IS_CHAIN_ATOM)
                    continue;
                struct xy_vector xy = moving_atom_offsets[atom_index++];
                add_collider(&list, board, (struct collider){
                    .center = {
                        v.x + (xy.x * rx - xy.y * ry),
                        v.y + (xy.x * ry + xy.y * rx),
                    },
                    .radius = atomRadius,
                });
            }
            while (chain_atom_cursor < list.number_of_chain_atom_colliders && list.chain_atom_colliders[chain_atom_cursor].movement_index == i) {
                list.chain_atom_colliders[chain_atom_cursor].v = v;
                list.chain_atom_colliders[chain_atom_cursor].r = (struct xy_vector){ rx, ry };
                if (board->area_growth_order == GROWTH_LINEAR)
                    mark_area_at_infinity(board, list.chain_atom_colliders[chain_atom_cursor]);
                chain_atom_cursor++;
            }
            list.cursor = list.length;
            list.bounding_box_up_to_cursor = list.bounding_box;
        }
        assert(chain_atom_cursor == list.number_of_chain_atom_colliders);
        if (chain_mode == EXTEND_CHAIN)
            resolve_chain_atom_collisions(&list);

        // printf("[");
        // for (size_t i = 0; i < list.length; ++i) {
        //     if (i > 0)
        //         printf(",");
        //     printf("[%f,%f,%f]", list.colliders[i].radius, list.colliders[i].center.x, list.colliders[i].center.y);
        // }
        // printf("],");
    }
    free(moving_atom_offsets);
    free(list.colliders);
    free(list.chain_atom_colliders);
    return list.collision;
}
