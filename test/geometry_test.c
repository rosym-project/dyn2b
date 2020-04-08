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

static struct ga_twist xda = {
    .target_body = &body_a, .reference_body = &body_b,
    .point = &point_a, .frame = &frame_a
};

static struct gc_twist xdc = {
    .angular_velocity = (struct vector3 [1]) { { 1.0, 2.0, 3.0 } },
    .linear_velocity = (struct vector3 [1]) { { 2.0, 3.0, 4.0 } }
};


START_TEST(test_gc_pose_compose)
{
    struct gc_pose x1 = {
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

    gc_pose_compose(&x1, &xc, &r);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(r.rotation->row[i].data[j], res_ang.row[i].data[j]);
        }
        ck_assert_flt_eq(r.translation->data[i], res_lin.data[i]);
    }


    struct gc_pose x2 = {
        .rotation = (struct matrix3x3 []) { {
                .row_x = {  cos(M_PI_4), -sin(M_PI_4), 0.0 },
                .row_y = {  0.0        ,  0.0        , 1.0 },
                .row_z = { -sin(M_PI_4), -cos(M_PI_4), 0.0 } } },
        .translation = (struct vector3 []) {
            (struct vector3) { 3.0, 2.0, 1.0 } }
    };

    struct matrix3x3 res_ang2 = {
        .row_x = { 0.0,  cos(M_PI_4), -sin(M_PI_4) },
        .row_y = { 1.0,  0.0        ,  0.0         },
        .row_z = { 0.0, -sin(M_PI_4), -cos(M_PI_4) } };
    struct vector3 res_lin2 = { 2.0, 5.0, 5.0 };

    gc_pose_compose(&x2, &xc, &r);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(r.rotation->row[i].data[j], res_ang2.row[i].data[j]);
        }
        ck_assert_flt_eq(r.translation->data[i], res_lin2.data[i]);
    }
}
END_TEST


START_TEST(test_ga_pose_compose)
{
    struct ga_pose x = {
        .target_body = &body_c, .target_frame = &frame_c,
        .reference_body = &body_b, .reference_frame = &frame_b };
    struct ga_pose r;

    ga_pose_compose(&x, &xa, &r);
    ck_assert_ptr_eq(r.target_body, &body_c);
    ck_assert_ptr_eq(r.target_frame, &frame_c);
    ck_assert_ptr_eq(r.reference_body, &body_a);
    ck_assert_ptr_eq(r.reference_frame, &frame_a);
}
END_TEST

START_TEST(test_gc_twist_tf_ref_to_tgt)
{
    struct gc_twist r = {
        .angular_velocity = (struct vector3 [1]) {},
        .linear_velocity = (struct vector3 [1]) {} };

    struct vector3 res_ang = { 2.0, 3.0, 1.0 };
    struct vector3 res_lin = { 3.0, 4.0, 2.0 };

    gc_twist_tf_ref_to_tgt(&xc, &xdc, &r);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(r.angular_velocity->data[i], res_ang.data[i]);
        ck_assert_flt_eq(r.linear_velocity->data[i], res_lin.data[i]);
    }
}
END_TEST


START_TEST(test_ga_twist_tf_ref_to_tgt)
{
    struct ga_twist r;

    ga_twist_tf_ref_to_tgt(&xa, &xda, &r);
    ck_assert_ptr_eq(r.target_body, &body_a);
    ck_assert_ptr_eq(r.reference_body, &body_b);
    ck_assert_ptr_eq(r.frame, &frame_b);
    ck_assert_ptr_eq(r.point, frame_b.origin);
}
END_TEST


START_TEST(test_gc_twist_accumulate)
{
    struct gc_twist r = {
        .angular_velocity = (struct vector3 [1]) {},
        .linear_velocity = (struct vector3 [1]) {} };

    struct vector3 res_ang = { 2.0, 4.0, 6.0 };
    struct vector3 res_lin = { 4.0, 6.0, 8.0 };

    gc_twist_accumulate(&xdc, &xdc, &r);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(r.angular_velocity->data[i], res_ang.data[i]);
        ck_assert_flt_eq(r.linear_velocity->data[i], res_lin.data[i]);
    }
}
END_TEST


START_TEST(test_ga_twist_accumulate)
{
    struct ga_twist xd = {
        .target_body = &body_c, .reference_body = &body_a,
        .point = &point_a, .frame = &frame_a };
    struct ga_twist r;

    ga_twist_accumulate(&xda, &xd, &r);
    ck_assert_ptr_eq(r.target_body, &body_c);
    ck_assert_ptr_eq(r.reference_body, &body_b);
    ck_assert_ptr_eq(r.frame, &frame_a);
    ck_assert_ptr_eq(r.point, frame_a.origin);
}
END_TEST


TCase *geometry_test()
{
    TCase *tc = tcase_create("Geometry");

    tcase_add_test(tc, test_gc_pose_compose);
    tcase_add_test(tc, test_ga_pose_compose);
    tcase_add_test(tc, test_gc_twist_tf_ref_to_tgt);
    tcase_add_test(tc, test_ga_twist_tf_ref_to_tgt);
    tcase_add_test(tc, test_gc_twist_accumulate);
    tcase_add_test(tc, test_ga_twist_accumulate);

    return tc;
}
