#include <dyn2b/functions/geometry.h>
#include <check.h>
#include <math.h>


#define ck_assert_flt_eq(X, Y) ck_assert_double_eq_tol(X, Y, 0.0001)


static struct point point_a;
static struct point point_b;
static struct point point_c;
static struct body body_a;
static struct body body_b;
static struct body body_c;
static struct frame frame_a = { .origin = &point_a };
static struct frame frame_b = { .origin = &point_b };
static struct frame frame_c = { .origin = &point_c };

static const struct ga_pose xa = {
    .target_body = &body_b, .target_frame = &frame_b,
    .reference_body = &body_a, .reference_frame = &frame_a
};

static struct gc_pose xc = {
    .rotation = (struct matrix3x3 [1]) { {
            .row_x = { 0.0, 1.0, 0.0 },
            .row_y = { 0.0, 0.0, 1.0 },
            .row_z = { 1.0, 0.0, 0.0 } } },
    .translation = (struct vector3 [1]) {
        { 1.0, 2.0, 3.0 } }
};


START_TEST(test_gc_pose_compose)
{
    struct gc_pose x = {
        .rotation = (struct matrix3x3 []) { {
                .row_x = { 0.0, 1.0, 0.0 },
                .row_y = { 0.0, 0.0, 1.0 },
                .row_z = { 1.0, 0.0, 0.0 } } },
        .translation = (struct vector3 []) {
            (struct vector3) { 1.0, 2.0, 3.0 } }
    };
    struct gc_pose r = {
        .rotation = (struct matrix3x3 [1]) {},
        .translation = (struct vector3 [1]) {}
    };

    struct matrix3x3 res_ang = {
        .row_x = { 0.0, 0.0, 1.0 },
        .row_y = { 1.0, 0.0, 0.0 },
        .row_z = { 0.0, 1.0, 0.0 } };
    struct vector3 res_lin = { 4.0, 3.0, 5.0 };

    gc_pose_compose(&xc, &x, &r);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(r.rotation->row[i].data[j], res_ang.row[i].data[j]);
        }
        ck_assert_flt_eq(r.translation->data[i], res_lin.data[i]);
    }
}
END_TEST


START_TEST(test_ga_pose_compose)
{
    struct ga_pose x = {
        .target_body = &body_c, .target_frame = &frame_c,
        .reference_body = &body_b, .reference_frame = &frame_b };
    struct ga_pose r;

    ga_pose_compose(&xa, &x, &r);
    ck_assert_ptr_eq(r.target_body, &body_c);
    ck_assert_ptr_eq(r.target_frame, &frame_c);
    ck_assert_ptr_eq(r.reference_body, &body_a);
    ck_assert_ptr_eq(r.reference_frame, &frame_a);
}
END_TEST


TCase *geometry_test()
{
    TCase *tc = tcase_create("Geometry");

    tcase_add_test(tc, test_gc_pose_compose);
    tcase_add_test(tc, test_ga_pose_compose);

    return tc;
}