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

struct gc_twist xdc = {
    .angular_velocity = (struct vector3 [1]) { 1.0, 2.0, 3.0 },
    .linear_velocity = (struct vector3 [1]) { 3.0, 4.0, 5.0 }
};

struct ga_twist xda = {
    .target_body = &body_a, .reference_body = &body_b,
    .point = &point_a, .frame = &frame_a
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

static struct ma_abi ma = {
    .body = &body_b, .point = &point_b, .frame = &frame_b
};

static struct mc_abi mc = {
    .zeroth_moment_of_mass = {
        .row_x = { 2.0, 0.0, 0.0 },
        .row_y = { 0.0, 3.0, 0.0 },
        .row_z = { 0.0, 0.0, 4.0 } },
    .first_moment_of_mass = {
        .row_x = { 4.0, 0.0, 0.0 },
        .row_y = { 0.0, 5.0, 0.0 },
        .row_z = { 0.0, 0.0, 6.0 } },
    .second_moment_of_mass = {
        .row_x = { 3.0, 4.0, 5.0 },
        .row_y = { 4.0, 6.0, 7.0 },
        .row_z = { 5.0, 7.0, 8.0 } }
};


START_TEST(test_mc_momentum_derive)
{
    struct mc_momentum p = {
        .angular_momentum = (struct vector3 [1]) { 24.0, 41.0, 41.0 },
        .linear_momentum = (struct vector3 [1]) { 4.0, 12.0,  8.0 } };
    struct mc_wrench r = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {} };

    struct vector3 res_ang = { -69.0, 27.0, 13.0 };
    struct vector3 res_lin = { -20.0,  4.0,  4.0 };

    mc_momentum_derive(&xdc, &p, &r);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(r.torque->data[i], res_ang.data[i]);
        ck_assert_flt_eq(r.force->data[i], res_lin.data[i]);
    }
}
END_TEST


START_TEST(test_ma_momentum_derive)
{
    struct ma_momentum p = {
        .body = &body_a,
        .point = &point_a,
        .frame = &frame_a };
    struct ma_wrench r;

    ma_momentum_derive(&xda, &p, &r);
    ck_assert_ptr_eq(r.body, &body_a);
    ck_assert_ptr_eq(r.point, &point_a);
    ck_assert_ptr_eq(r.frame, &frame_a);
}
END_TEST


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


START_TEST(test_mc_wrench_invert)
{
    struct mc_wrench f = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {} };
    struct vector3 res_ang[2] = { { -1.0, -2.0, -3.0 }, { -2.0, - 4.0, - 6.0 } };
    struct vector3 res_lin[2] = { { -2.0, -3.0, -4.0 }, { -8.0, -10.0, -12.0 } };

    mc_wrench_invert(&fc, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin[j].data[i]);
        }
    }
}
END_TEST


START_TEST(test_ma_wrench_invert)
{
    struct ma_wrench f;

    ma_wrench_invert(&fa, &f);
    ck_assert_ptr_eq(f.body, &body_b);
    ck_assert_ptr_eq(f.point, frame_b.origin);
    ck_assert_ptr_eq(f.frame, &frame_b);
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


START_TEST(test_mc_rbi_map_twist_to_momentum)
{
    struct mc_rbi m = {
        .zeroth_moment_of_mass = 2.0,
        .first_moment_of_mass = { 4.0, 6.0, 8.0 },
        .second_moment_of_mass = {
            .row_x = { 3.0, 4.0, 5.0 },
            .row_y = { 4.0, 6.0, 7.0 },
            .row_z = { 5.0, 7.0, 8.0 } } };
    struct mc_momentum r = {
        .angular_momentum = (struct vector3 [1]) {},
        .linear_momentum = (struct vector3 [1]) {} };

    struct vector3 res_ang = { 24.0, 41.0, 41.0 };
    struct vector3 res_lin = {  4.0, 12.0,  8.0 };

    mc_rbi_map_twist_to_momentum(&m, &xdc, &r);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(r.angular_momentum[0].data[i], res_ang.data[i]);
        ck_assert_flt_eq(r.linear_momentum[0].data[i], res_lin.data[i]);
    }
}
END_TEST


START_TEST(test_ma_rbi_map_twist_to_momentum)
{
    struct ma_rbi m = {
        .body = &body_a,
        .point = &point_a,
        .frame = &frame_a
    };
    struct ma_momentum r;

    ma_rbi_map_twist_to_momentum(&m, &xda, &r);
    ck_assert_ptr_eq(r.body, &body_a);
    ck_assert_ptr_eq(r.point, &point_a);
    ck_assert_ptr_eq(r.frame, &frame_a);
}
END_TEST


START_TEST(test_mc_rbi_to_abi)
{
    struct mc_rbi m = {
        .zeroth_moment_of_mass = 2.0,
        .first_moment_of_mass = { 4.0, 6.0, 8.0 },
        .second_moment_of_mass = {
            .row_x = { 3.0, 4.0, 5.0 },
            .row_y = { 4.0, 6.0, 7.0 },
            .row_z = { 5.0, 7.0, 8.0 } } };
    struct mc_abi r;

    struct matrix3x3 res_m2 = {
            .row_x = { 3.0, 4.0, 5.0 },
            .row_y = { 4.0, 6.0, 7.0 },
            .row_z = { 5.0, 7.0, 8.0 } };
    struct matrix3x3 res_m1 = {
        .row_x = {  0.0, -8.0,  6.0 },
        .row_y = {  8.0,  0.0, -4.0 },
        .row_z = { -6.0,  4.0,  0.0 } };
    struct matrix3x3 res_m0 = {
        .row_x = { 2.0, 0.0, 0.0 },
        .row_y = { 0.0, 2.0, 0.0 },
        .row_z = { 0.0, 0.0, 2.0 } };

    mc_rbi_to_abi(&m, &r);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(r.second_moment_of_mass.row[i].data[j], res_m2.row[i].data[j]);
            ck_assert_flt_eq(r.first_moment_of_mass.row[i].data[j], res_m1.row[i].data[j]);
            ck_assert_flt_eq(r.zeroth_moment_of_mass.row[i].data[j], res_m0.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_ma_rbi_to_abi)
{
    struct ma_rbi m = {
        .body = &body_a,
        .point = &point_a,
        .frame = &frame_a
    };
    struct ma_abi r;

    ma_rbi_to_abi(&m, &r);
    ck_assert_ptr_eq(r.body, &body_a);
    ck_assert_ptr_eq(r.point, &point_a);
    ck_assert_ptr_eq(r.frame, &frame_a);
}
END_TEST


START_TEST(test_mc_abi_tf_tgt_to_ref)
{
    struct mc_abi m;

    struct matrix3x3 res_m = {
        .row_x = { 4.0, 0.0, 0.0 },
        .row_y = { 0.0, 2.0, 0.0 },
        .row_z = { 0.0, 0.0, 3.0 } };
    struct matrix3x3 res_h = {
        .row_x = {   6.0, -6.0,  6.0 },
        .row_y = {  12.0,  4.0, -3.0 },
        .row_z = { - 8.0,  2.0,  5.0 } };
    struct matrix3x3 res_i = {
        .row_x = {  38.0,   5.0, - 1.0 },
        .row_y = {   5.0,  42.0, -21.0 },
        .row_z = { - 1.0, -21.0,  24.0 } };

    mc_abi_tf_tgt_to_ref(&xc, &mc, &m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(m.zeroth_moment_of_mass.row[i].data[j], res_m.row[i].data[j]);
            ck_assert_flt_eq(m.first_moment_of_mass.row[i].data[j], res_h.row[i].data[j]);
            ck_assert_flt_eq(m.second_moment_of_mass.row[i].data[j], res_i.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_ma_abi_tf_tgt_to_ref)
{
    struct ma_abi m;

    ma_abi_tf_tgt_to_ref(&xa, &ma, &m);
    ck_assert_ptr_eq(m.point, frame_a.origin);
    ck_assert_ptr_eq(m.body, &body_b);
    ck_assert_ptr_eq(m.frame, &frame_a);
}
END_TEST


START_TEST(test_mc_abi_add)
{
    struct mc_abi r;

    struct matrix3x3 res_m2 = {
        .row_x = {  6.0,  8.0, 10.0 },
        .row_y = {  8.0, 12.0, 14.0 },
        .row_z = { 10.0, 14.0, 16.0 } };
    struct matrix3x3 res_m1 = {
        .row_x = { 8.0,  0.0, 0.0 },
        .row_y = { 0.0, 10.0, 0.0 },
        .row_z = { 0.0,  0.0, 12.0 } };
    struct matrix3x3 res_m0 = {
        .row_x = { 4.0, 0.0, 0.0 },
        .row_y = { 0.0, 6.0, 0.0 },
        .row_z = { 0.0, 0.0, 8.0 } };

    mc_abi_add(&mc, &mc, &r);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(r.second_moment_of_mass.row[i].data[j], res_m2.row[i].data[j]);
            ck_assert_flt_eq(r.first_moment_of_mass.row[i].data[j], res_m1.row[i].data[j]);
            ck_assert_flt_eq(r.zeroth_moment_of_mass.row[i].data[j], res_m0.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_ma_abi_add)
{
    struct ma_abi m;

    ma_abi_add(&ma, &ma, &m);
    ck_assert_ptr_eq(m.body, &body_b);
    ck_assert_ptr_eq(m.point, &point_b);
    ck_assert_ptr_eq(m.frame, &frame_b);
}
END_TEST


TCase *mechanics_test()
{
    TCase *tc = tcase_create("Mechanics");

    tcase_add_test(tc, test_mc_momentum_derive);
    tcase_add_test(tc, test_ma_momentum_derive);
    tcase_add_test(tc, test_mc_wrench_tf_tgt_to_ref);
    tcase_add_test(tc, test_ma_wrench_tf_tgt_to_ref);
    tcase_add_test(tc, test_mc_wrench_invert);
    tcase_add_test(tc, test_ma_wrench_invert);
    tcase_add_test(tc, test_mc_wrench_add);
    tcase_add_test(tc, test_ma_wrench_add);
    tcase_add_test(tc, test_mc_wrench_sub);
    tcase_add_test(tc, test_ma_wrench_sub);
    tcase_add_test(tc, test_mc_rbi_map_twist_to_momentum);
    tcase_add_test(tc, test_ma_rbi_map_twist_to_momentum);
    tcase_add_test(tc, test_mc_rbi_to_abi);
    tcase_add_test(tc, test_ma_rbi_to_abi);
    tcase_add_test(tc, test_mc_abi_tf_tgt_to_ref);
    tcase_add_test(tc, test_ma_abi_tf_tgt_to_ref);
    tcase_add_test(tc, test_mc_abi_add);
    tcase_add_test(tc, test_ma_abi_add);

    return tc;
}
