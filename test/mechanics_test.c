#include <dyn2b/functions/mechanics.h>
#include <check.h>
#include <math.h>


#define ck_assert_flt_eq(X, Y) ck_assert_double_eq_tol(X, Y, 0.0001)


static struct point point_a;
static struct point point_b;
static struct body body_a;
static struct body body_b;
static struct frame frame_a = { .origin = &point_a };
static struct frame frame_b = { .origin = &point_b };

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

static const struct ma_wrench fa = {
    .body = &body_b,
    .point = &point_b,
    .frame = &frame_b
};

static const struct mc_wrench fc = {
    .torque = (struct vector3 [2]) {
        { 1.0, 2.0, 3.0 },
        { 2.0, 4.0, 6.0 } },
    .force = (struct vector3 [2]) {
        { 2.0,  3.0,  4.0 },
        { 8.0, 10.0, 12.0 } }
};


START_TEST(test_mc_wrench_tf_tgt_to_ref)
{
    struct mc_wrench f = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {} };
    struct vector3 res_ang[2] = { { 3.0, 10.0, -4.0 }, {  2.0, 28.0, -12.0 } };
    struct vector3 res_lin[2] = { { 4.0,  2.0,  3.0 }, { 12.0,  8.0,  10.0 } };

    mc_wrench_tf_tgt_to_ref(&xc, &fc, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin[j].data[i]);
        }
    }
}
END_TEST


START_TEST(test_ma_wrench_tf_tgt_to_ref)
{
    struct ma_wrench f;

    ma_wrench_tf_tgt_to_ref(&xa, &fa, &f);
    ck_assert_ptr_eq(f.body, &body_b);
    ck_assert_ptr_eq(f.point, frame_a.origin);
    ck_assert_ptr_eq(f.frame, &frame_a);
}
END_TEST


START_TEST(test_mc_wrench_add)
{
    struct mc_wrench f = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {} };
    struct vector3 res_ang[2] = { { 2.0, 4.0, 6.0 }, {  4.0,  8.0, 12.0 } };
    struct vector3 res_lin[2] = { { 4.0, 6.0, 8.0 }, { 16.0, 20.0, 24.0 } };

    mc_wrench_add(&fc, &fc, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin[j].data[i]);
        }
    }
}
END_TEST


START_TEST(test_ma_wrench_add)
{
    struct ma_wrench f;

    ma_wrench_add(&fa, &fa, &f);
    ck_assert_ptr_eq(f.body, &body_b);
    ck_assert_ptr_eq(f.point, frame_b.origin);
    ck_assert_ptr_eq(f.frame, &frame_b);
}
END_TEST


START_TEST(test_mc_wrench_sub)
{
    struct mc_wrench f = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {} };
    struct vector3 res_ang[2] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
    struct vector3 res_lin[2] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };

    mc_wrench_sub(&fc, &fc, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin[j].data[i]);
        }
    }
}
END_TEST


START_TEST(test_ma_wrench_sub)
{
    struct ma_wrench f;

    ma_wrench_sub(&fa, &fa, &f);
    ck_assert_ptr_eq(f.body, &body_b);
    ck_assert_ptr_eq(f.point, frame_b.origin);
    ck_assert_ptr_eq(f.frame, &frame_b);
}
END_TEST


TCase *mechanics_test()
{
    TCase *tc = tcase_create("Mechanics");

    tcase_add_test(tc, test_mc_wrench_tf_tgt_to_ref);
    tcase_add_test(tc, test_ma_wrench_tf_tgt_to_ref);
    tcase_add_test(tc, test_mc_wrench_add);
    tcase_add_test(tc, test_ma_wrench_add);
    tcase_add_test(tc, test_mc_wrench_sub);
    tcase_add_test(tc, test_ma_wrench_sub);

    return tc;
}