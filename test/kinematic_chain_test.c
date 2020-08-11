#include <dyn2b/functions/kinematic_chain.h>
#include <check.h>
#include <math.h>


#ifdef ck_assert_double_eq_tol
#  define ck_assert_flt_eq(X, Y) ck_assert_double_eq_tol(X, Y, 0.0001)
#else
#  define ck_assert_flt_eq(X, Y) do { \
     double _dist = fabs((double)(X) - (double)(Y)); \
     ck_assert_msg(_dist < (0.0001), "Assertion '%s' failed: %s == %f, %s == %f", #X" == "#Y, #X, (X), #Y, (Y)); \
   } while (0)
#endif


static struct body body_a;
static struct body body_b;
static struct body body_c;
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

static struct ga_twist xda = {
    .target_body = &body_b, .reference_body = &body_c,
    .point = &point_b, .frame = &frame_b
};

static struct gc_twist xdc = {
    .angular_velocity = (struct vector3 [1]) { { 1.0, 2.0, 3.0 } },
    .linear_velocity = (struct vector3 [1]) { { 2.0, 3.0, 4.0 } }
};

static struct ma_abi ma = {
    .body = &body_b,
    .point = &point_b,
    .frame = &frame_b
};

static struct mc_abi mc = {
    .zeroth_moment_of_mass = {
        .row_x = { 4.0, 5.0, 6.0 },
        .row_y = { 5.0, 5.0, 6.0 },
        .row_z = { 6.0, 6.0, 6.0 }
    },
    .first_moment_of_mass = {
        .row_x = { 2.0, 3.0, 4.0 },
        .row_y = { 3.0, 3.0, 4.0 },
        .row_z = { 4.0, 4.0, 4.0 }
    },
    .second_moment_of_mass = {
        .row_x = { 1.0, 2.0, 3.0 },
        .row_y = { 2.0, 2.0, 3.0 },
        .row_z = { 3.0, 3.0, 3.0 }
    }
};

static struct ma_wrench fa = {
    .body = &body_b,
    .point = &point_b,
    .frame = &frame_b
};

static struct mc_wrench fc = {
    .torque = (struct vector3 [2]) {
        { 1.0, 2.0, 3.0 },
        { 2.0, 4.0, 6.0 } },
    .force = (struct vector3 [2]) {
        { 2.0,  3.0,  4.0 },
        { 8.0, 10.0, 12.0 } }
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


START_TEST(test_kca_fvk)
{
    struct ga_twist xd;

    kca_fvk(&ja, &xd);
    ck_assert_ptr_eq(xd.target_body, &body_b);
    ck_assert_ptr_eq(xd.reference_body, &body_a);
    ck_assert_ptr_eq(xd.point, &point_b);
    ck_assert_ptr_eq(xd.frame, &frame_b);
}
END_TEST


START_TEST(test_kca_fak)
{
    struct ga_acc_twist xdd;

    kca_fak(&ja, &xdd);
    ck_assert_ptr_eq(xdd.target_body, &body_b);
    ck_assert_ptr_eq(xdd.reference_body, &body_a);
    ck_assert_ptr_eq(xdd.frame, &frame_b);
    ck_assert_ptr_eq(xdd.point, &point_b);
}
END_TEST


START_TEST(test_kca_inertial_acceleration)
{
    struct ga_acc_twist xdd;

    kca_inertial_acceleration(&ja, &xda, &xdd);
    ck_assert_ptr_eq(xdd.target_body, &body_b);
    ck_assert_ptr_eq(xdd.reference_body, &body_a);
    ck_assert_ptr_eq(xdd.frame, &frame_b);
    ck_assert_ptr_eq(xdd.point, &point_b);
}
END_TEST


START_TEST(test_kca_ifk)
{
    kca_ifk(&ja, &fa);
}
END_TEST


START_TEST(test_kca_ffd)
{
    struct ma_wrench r;

    kca_ffd(&ja, &ma, &r);
    ck_assert_ptr_eq(r.body, &body_a);
    ck_assert_ptr_eq(r.frame, &frame_b);
    ck_assert_ptr_eq(r.point, &point_b);
}
END_TEST


START_TEST(test_kca_project_inertia)
{
    struct ma_abi m;

    kca_project_inertia(&ja, &ma, &m);
    ck_assert_ptr_eq(m.body, &body_a);
    ck_assert_ptr_eq(m.point, &point_b);
    ck_assert_ptr_eq(m.frame, &frame_b);
}
END_TEST


START_TEST(test_kca_project_wrench)
{
    struct ma_wrench f;

    kca_project_wrench(&ja, &ma, &fa, &f);
    ck_assert_ptr_eq(f.body, &body_a);
    ck_assert_ptr_eq(f.point, &point_b);
    ck_assert_ptr_eq(f.frame, &frame_b);
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


START_TEST(test_rev_fvk)
{
    struct kcc_joint joint = { .type = JOINT_TYPE_REVOLUTE };
    joint_velocity qd = { 2.0 };
    struct gc_twist xd = {
        .angular_velocity = (struct vector3 [1]) {},
        .linear_velocity = (struct vector3 [1]) {}
    };

    struct vector3 res_ang_x = { 2.0, 0.0, 0.0 };
    struct vector3 res_lin_x = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].fvk(&joint, &qd, &xd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xd.angular_velocity->data[i], res_ang_x.data[i]);
        ck_assert_flt_eq(xd.linear_velocity->data[i], res_lin_x.data[i]);
    }


    struct vector3 res_ang_y = { 0.0, 2.0, 0.0 };
    struct vector3 res_lin_y = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].fvk(&joint, &qd, &xd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xd.angular_velocity->data[i], res_ang_y.data[i]);
        ck_assert_flt_eq(xd.linear_velocity->data[i], res_lin_y.data[i]);
    }


    struct vector3 res_ang_z = { 0.0, 0.0, 2.0 };
    struct vector3 res_lin_z = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].fvk(&joint, &qd, &xd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xd.angular_velocity->data[i], res_ang_z.data[i]);
        ck_assert_flt_eq(xd.linear_velocity->data[i], res_lin_z.data[i]);
    }
}
END_TEST


START_TEST(test_rev_fak)
{
    struct kcc_joint joint = { .type = JOINT_TYPE_REVOLUTE };
    joint_acceleration qdd = { 2.0 };
    struct gc_acc_twist xdd = {
        .angular_acceleration = (struct vector3 [1]) {},
        .linear_acceleration = (struct vector3 [1]) {}
    };

    struct vector3 res_ang_x = { 2.0, 0.0, 0.0 };
    struct vector3 res_lin_x = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].fak(&joint, &qdd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_x.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_x.data[i]);
    }


    struct vector3 res_ang_y = { 0.0, 2.0, 0.0 };
    struct vector3 res_lin_y = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].fak(&joint, &qdd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_y.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_y.data[i]);
    }


    struct vector3 res_ang_z = { 0.0, 0.0, 2.0 };
    struct vector3 res_lin_z = { 0.0, 0.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].fak(&joint, &qdd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_z.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_z.data[i]);
    }
}
END_TEST


START_TEST(test_rev_inertial_acceleration)
{
    struct kcc_joint joint = { .type = JOINT_TYPE_REVOLUTE };
    joint_position qd = { 2.0 };
    struct gc_twist xd = {
        .angular_velocity = (struct vector3 [1]) { 2.0, 3.0, 4.0 },
        .linear_velocity = (struct vector3 [1]) { 3.0, 4.0, 5.0 }
    };
    struct gc_acc_twist xdd = {
        .angular_acceleration = (struct vector3 [1]) {},
        .linear_acceleration = (struct vector3 [1]) {}
    };


    struct vector3 res_ang_x = { 0.0,  8.0, -6.0 };
    struct vector3 res_lin_x = { 0.0, 10.0, -8.0 };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].inertial_acceleration(&joint, &xd, &qd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_x.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_x.data[i]);
    }


    struct vector3 res_ang_y = { - 8.0, 0.0, 4.0 };
    struct vector3 res_lin_y = { -10.0, 0.0, 6.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].inertial_acceleration(&joint, &xd, &qd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_y.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_y.data[i]);
    }


    struct vector3 res_ang_z = { 6.0, -4.0, 0.0 };
    struct vector3 res_lin_z = { 8.0, -6.0, 0.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].inertial_acceleration(&joint, &xd, &qd, &xdd);
    for (int i = 0; i < 3; i++) {
        ck_assert_flt_eq(xdd.angular_acceleration->data[i], res_ang_z.data[i]);
        ck_assert_flt_eq(xdd.linear_acceleration->data[i], res_lin_z.data[i]);
    }
}
END_TEST


START_TEST(test_rev_ifk)
{
    struct kcc_joint joint = { .type = JOINT_TYPE_REVOLUTE };
    joint_torque tau[2] = { 0.0, 0.0 };

    joint_torque res_x[2] = { 1.0, 2.0 };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].ifk(&joint, &fc, tau, 2);
    for (int i = 0; i < 2; i++) ck_assert_flt_eq(tau[i], res_x[i]);


    joint_torque res_y[2] = { 2.0, 4.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].ifk(&joint, &fc, tau, 2);
    for (int i = 0; i < 2; i++) ck_assert_flt_eq(tau[i], res_y[i]);


    joint_torque res_z[2] = { 3.0, 6.0 };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].ifk(&joint, &fc, tau, 2);
    for (int i = 0; i < 2; i++) ck_assert_flt_eq(tau[i], res_z[i]);
}
END_TEST


START_TEST(test_rev_ffd)
{
    struct kcc_joint joint = {
        .type = JOINT_TYPE_REVOLUTE,
        .revolute_joint.inertia = (double [1]) { 3.0 }
    };
    joint_torque tau[2] = { 1.0, 2.0 };
    struct mc_wrench f = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {}
    };


    struct vector3 res_ang_x[2] = { { 0.25, 0.5 , 0.75 }, { 0.5, 1.0, 1.5 } };
    struct vector3 res_lin_x[2] = { { 0.5 , 0.75, 1.0  }, { 1.0, 1.5, 2.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].ffd(&joint, &mc, tau, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang_x[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin_x[j].data[i]);
        }
    }


    struct vector3 res_ang_y[2] = { { 0.4, 0.4, 0.6 }, { 0.8, 0.8, 1.2 } };
    struct vector3 res_lin_y[2] = { { 0.6, 0.6, 0.8 }, { 1.2, 1.2, 1.6 } };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].ffd(&joint, &mc, tau, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang_y[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin_y[j].data[i]);
        }
    }


    struct vector3 res_ang_z[2] = {
            { 0.5, 0.5, 0.5 },
            { 1.0, 1.0, 1.0 } };
    struct vector3 res_lin_z[2] = {
            { 2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0 },
            { 4.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].ffd(&joint, &mc, tau, &f, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(f.torque[j].data[i], res_ang_z[j].data[i]);
            ck_assert_flt_eq(f.force[j].data[i], res_lin_z[j].data[i]);
        }
    }
}
END_TEST


START_TEST(test_rev_project_inertia)
{
    struct kcc_joint joint = {
        .type = JOINT_TYPE_REVOLUTE,
        .revolute_joint.inertia = (double [1]) { 3.0 }
    };
    struct mc_abi m;

    struct matrix3x3 res_m2_x = {
        .row_x = { 0.75, 1.5, 2.25 },
        .row_y = { 1.5 , 1.0, 1.5  },
        .row_z = { 2.25, 1.5, 0.75 } };
    struct matrix3x3 res_m1_x = {
        .row_x = { 1.5, 2.25, 3.0 },
        .row_y = { 2.0, 1.5 , 2.0 },
        .row_z = { 2.5, 1.75, 1.0 } };
    struct matrix3x3 res_m0_x = {
        .row_x = { 3.0, 3.5 , 4.0 },
        .row_y = { 3.5, 2.75, 3.0 },
        .row_z = { 4.0, 3.0 , 2.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_inertia(&joint, &mc, &m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(m.second_moment_of_mass.row[i].data[j], res_m2_x.row[i].data[j]);
            ck_assert_flt_eq(m.first_moment_of_mass.row[i].data[j], res_m1_x.row[i].data[j]);
            ck_assert_flt_eq(m.zeroth_moment_of_mass.row[i].data[j], res_m0_x.row[i].data[j]);
        }
    }

    struct matrix3x3 res_m2_y = {
        .row_x = { 0.2, 1.2, 1.8 },
        .row_y = { 1.2, 1.2, 1.8 },
        .row_z = { 1.8, 1.8, 1.2 } };
    struct matrix3x3 res_m1_y = {
        .row_x = { 0.8, 1.8, 2.4 },
        .row_y = { 1.8, 1.8, 2.4 },
        .row_z = { 2.2, 2.2, 1.6 } };
    struct matrix3x3 res_m0_y = {
        .row_x = { 2.2, 3.2, 3.6 },
        .row_y = { 3.2, 3.2, 3.6 },
        .row_z = { 3.6, 3.6, 2.8 } };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_inertia(&joint, &mc, &m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(m.second_moment_of_mass.row[i].data[j], res_m2_y.row[i].data[j]);
            ck_assert_flt_eq(m.first_moment_of_mass.row[i].data[j], res_m1_y.row[i].data[j]);
            ck_assert_flt_eq(m.zeroth_moment_of_mass.row[i].data[j], res_m0_y.row[i].data[j]);
        }
    }


    struct matrix3x3 res_m2_z = {
        .row_x = { -0.5, 0.5, 1.5 },
        .row_y = {  0.5, 0.5, 1.5 },
        .row_z = {  1.5, 1.5, 1.5 } };
    struct matrix3x3 res_m1_z = {
        .row_x = { 0.0, 1.0, 2.0 },
        .row_y = { 1.0, 1.0, 2.0 },
        .row_z = { 2.0, 2.0, 2.0 } };
    struct matrix3x3 res_m0_z = {
        .row_x = {  4.0 / 3.0,  7.0 / 3.0, 10.0 / 3.0 },
        .row_y = {  7.0 / 3.0,  7.0 / 3.0, 10.0 / 3.0 },
        .row_z = { 10.0 / 3.0, 10.0 / 3.0, 10.0 / 3.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_inertia(&joint, &mc, &m);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ck_assert_flt_eq(m.second_moment_of_mass.row[i].data[j], res_m2_z.row[i].data[j]);
            ck_assert_flt_eq(m.first_moment_of_mass.row[i].data[j], res_m1_z.row[i].data[j]);
            ck_assert_flt_eq(m.zeroth_moment_of_mass.row[i].data[j], res_m0_z.row[i].data[j]);
        }
    }
}
END_TEST


START_TEST(test_rev_project_wrench)
{
    struct kcc_joint joint = {
        .type = JOINT_TYPE_REVOLUTE,
        .revolute_joint.inertia = (double [1]) { 3.0 }
    };
    struct mc_wrench r = {
        .torque = (struct vector3 [2]) {},
        .force = (struct vector3 [2]) {}
    };


    struct vector3 res_ang_x[2] = { {  0.75, 1.5 , 2.25 }, {  1.5, 3.0,  4.5 } };
    struct vector3 res_lin_x[2] = { {  1.5 , 2.25, 3.0  }, {  7.0, 8.5, 10.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_X;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_wrench(&joint, &mc, &fc, &r, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(r.torque[j].data[i], res_ang_x[j].data[i]);
            ck_assert_flt_eq(r.force[j].data[i], res_lin_x[j].data[i]);
        }
    }

    struct vector3 res_ang_y[2] = { {  0.2 , 1.2 , 1.8  }, {  0.4, 2.4,  3.6 } };
    struct vector3 res_lin_y[2] = { {  0.8 , 1.8 , 2.4  }, {  5.6, 7.6,  8.8 } };

    joint.revolute_joint.axis = JOINT_AXIS_Y;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_wrench(&joint, &mc, &fc, &r, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(r.torque[j].data[i], res_ang_y[j].data[i]);
            ck_assert_flt_eq(r.force[j].data[i], res_lin_y[j].data[i]);
        }
    }


    struct vector3 res_ang_z[2] = { { -0.5 , 0.5 , 1.5  }, { -1.0, 1.0,  3.0 } };
    struct vector3 res_lin_z[2] = { {  0.0 , 1.0 , 2.0  }, {  4.0, 6.0,  8.0 } };

    joint.revolute_joint.axis = JOINT_AXIS_Z;
    kcc_joint[JOINT_TYPE_REVOLUTE].project_wrench(&joint, &mc, &fc, &r, 2);
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            ck_assert_flt_eq(r.torque[j].data[i], res_ang_z[j].data[i]);
            ck_assert_flt_eq(r.force[j].data[i], res_lin_z[j].data[i]);
        }
    }
}
END_TEST


TCase *kinematic_chain_test()
{
    TCase *tc = tcase_create("KinematicChain");

    tcase_add_test(tc, test_kca_fpk);
    tcase_add_test(tc, test_kca_fvk);
    tcase_add_test(tc, test_kca_fak);
    tcase_add_test(tc, test_kca_inertial_acceleration);
    tcase_add_test(tc, test_kca_ifk);
    tcase_add_test(tc, test_kca_ffd);
    tcase_add_test(tc, test_kca_project_inertia);
    tcase_add_test(tc, test_kca_project_wrench);
    tcase_add_test(tc, test_rev_fpk);
    tcase_add_test(tc, test_rev_fvk);
    tcase_add_test(tc, test_rev_fak);
    tcase_add_test(tc, test_rev_inertial_acceleration);
    tcase_add_test(tc, test_rev_ifk);
    tcase_add_test(tc, test_rev_ffd);
    tcase_add_test(tc, test_rev_project_inertia);
    tcase_add_test(tc, test_rev_project_wrench);

    return tc;
}
