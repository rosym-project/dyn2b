#include <dyn2b/functions/kinematic_chain.h>
#include <check.h>
#include <math.h>


#define ck_assert_flt_eq(X, Y) ck_assert_double_eq_tol(X, Y, 0.0001)


static struct body body_a;
static struct body body_b;
static struct point point_a;
static struct point point_b;
static struct frame frame_a = { .origin = &point_a };
static struct frame frame_b = { .origin = &point_b };

static struct kca_joint ja = {
    .target_body = &body_b,
    .target_frame = &frame_b,
    .reference_body = &body_a,
    .reference_frame = &frame_a
};


START_TEST(test_kca_fpk)
{
    struct ga_pose x;

    kca_fpk(&ja, &x);
    ck_assert_ptr_eq(x.target_body, &body_b);
    ck_assert_ptr_eq(x.reference_body, &body_a);
    ck_assert_ptr_eq(x.target_frame, &frame_b);
    ck_assert_ptr_eq(x.reference_frame, &frame_a);
}
END_TEST


START_TEST(test_rev_fpk)
{
    struct kcc_joint joint = { .type = JOINT_TYPE_REVOLUTE };
    joint_position q = { M_PI_2 };
    struct gc_pose x = {
        .rotation = (struct matrix3x3 [1]) { {
            .row_x = { 2.0, 2.0, 2.0 },
            .row_y = { 2.0, 2.0, 2.0 },
            .row_z = { 2.0, 2.0, 2.0 } } },
        .translation = (struct vector3 [1]) { { 2.0, 2.0, 2.0 } }
    };

    struct matrix3x3 res_ang_x = {
        .row_x = { 1.0,  0.0, 0.0 },
        .row_y = { 0.0,  0.0, 1.0 },
        .row_z = { 0.0, -1.0, 0.0 } };
    struct vector3 res_lin_x = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].fpk(&joint, &q, &x);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(x.rotation->row[i].data[j], res_ang_x.row[i].data[j]);
        }
        ck_assert_flt_eq(x.translation->data[i], res_lin_x.data[i]);
    }


    struct matrix3x3 res_ang_y = {
        .row_x = { 0.0, 0.0, -1.0 },
        .row_y = { 0.0, 1.0,  0.0 },
        .row_z = { 1.0, 0.0,  0.0 } };
    struct vector3 res_lin_y = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].fpk(&joint, &q, &x);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(x.rotation->row[i].data[j], res_ang_y.row[i].data[j]);
        }
        ck_assert_flt_eq(x.translation->data[i], res_lin_y.data[i]);
    }


    struct matrix3x3 res_ang_z = {
        .row_x = {  0.0, 1.0, 0.0 },
        .row_y = { -1.0, 0.0, 0.0 },
        .row_z = {  0.0, 0.0, 1.0 } };
    struct vector3 res_lin_z = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].fpk(&joint, &q, &x);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(x.rotation->row[i].data[j], res_ang_z.row[i].data[j]);
        }
        ck_assert_flt_eq(x.translation->data[i], res_lin_z.data[i]);
    }
}
END_TEST



TCase *kinematic_chain_test()
{
    TCase *tc = tcase_create("KinematicChain");

    tcase_add_test(tc, test_kca_fpk);
    tcase_add_test(tc, test_rev_fpk);

    return tc;
}