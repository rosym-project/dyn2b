#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/functions/geometry.h>
#include <math.h>
#include <string.h>
#include <assert.h>


void kca_fpk(
        const struct kca_joint *joint,
        struct ga_pose *x)
{
    assert(joint);
    assert(x);

    x->target_body = joint->target_body;
    x->reference_body = joint->reference_body;
    x->target_frame = joint->target_frame;
    x->reference_frame = joint->reference_frame;
}


void kca_fvk(
        const struct kca_joint *joint,
        struct ga_twist *xd)
{
    assert(joint);
    assert(xd);
    assert(joint->target_frame);

    xd->target_body = joint->target_body;
    xd->reference_body = joint->reference_body;
    xd->point = joint->target_frame->origin;
    xd->frame = joint->target_frame;
}


void kca_fak(
        const struct kca_joint *joint,
        struct ga_acc_twist *xdd)
{
    assert(joint);
    assert(xdd);

    xdd->target_body = joint->target_body;
    xdd->reference_body = joint->reference_body;
    xdd->point = joint->target_frame->origin;
    xdd->frame = joint->target_frame;
}


void kca_inertial_acceleration(
        const struct kca_joint *joint,
        const struct ga_twist *xd,
        struct ga_acc_twist *xdd)
{
    assert(joint);
    assert(xd);
    assert(xdd);
    assert(xd->frame);
    assert(xd->frame->origin);
    assert(xd->point == xd->frame->origin);         // screw twist
    assert(xd->target_body == joint->target_body);
    assert(xd->frame == joint->target_frame);

    xdd->target_body = joint->target_body;
    xdd->reference_body = joint->reference_body;
    xdd->frame = joint->target_frame;
    xdd->point = joint->target_frame->origin;
}


static void rev_fpk(
        const struct kcc_joint *joint,
        const joint_position *q,
        struct gc_pose *x)
{
    assert(joint);
    assert(q);
    assert(x);
    assert(x->rotation);
    assert(x->translation);

    double cq = cos(q[0]);
    double sq = sin(q[0]);

    // initialize all rotation matrix and position vector to 0.0
    memset(x->rotation, 0, sizeof(*x->rotation));
    memset(x->translation, 0, sizeof(*x->translation));

    // Note that the rotation for spatial transforms is inverted when compared
    // to homogeneous transforms!
    if (joint->revolute_joint.axis == JOINT_AXIS_X) {
        // |1  0  0 |
        // |0  cq sq|
        // |0 -sq cq|

        x->rotation->row_x.x = 1.0;
        x->rotation->row_y.y =  cq;
        x->rotation->row_y.z =  sq;
        x->rotation->row_z.y = -sq;
        x->rotation->row_z.z =  cq;
    } else if (joint->revolute_joint.axis == JOINT_AXIS_Y) {
        // | cq 0 -sq|
        // | 0  1  0 |
        // | sq 0  cq|

        x->rotation->row_y.y = 1.0;
        x->rotation->row_x.x =  cq;
        x->rotation->row_x.z = -sq;
        x->rotation->row_z.x =  sq;
        x->rotation->row_z.z =  cq;
    } else if (joint->revolute_joint.axis == JOINT_AXIS_Z) {
        // | cq sq 0|
        // |-sq cq 0|
        // | 0  0  1|

        x->rotation->row_z.z = 1.0;
        x->rotation->row_x.x =  cq;
        x->rotation->row_x.y =  sq;
        x->rotation->row_y.x = -sq;
        x->rotation->row_y.y =  cq;
    } else {
        assert(0);
    }
}


static void rev_fvk(
        const struct kcc_joint *joint,
        const joint_velocity *qd,
        struct gc_twist *xd)
{
    assert(joint);
    assert(qd);
    assert(xd);
    assert(xd->angular_velocity);

    memset(xd->linear_velocity, 0, sizeof(*xd->linear_velocity));

    enum joint_axis axis = joint->revolute_joint.axis;
    xd->angular_velocity->x = (axis == JOINT_AXIS_X) ? qd[0] : 0.0;
    xd->angular_velocity->y = (axis == JOINT_AXIS_Y) ? qd[0] : 0.0;
    xd->angular_velocity->z = (axis == JOINT_AXIS_Z) ? qd[0] : 0.0;
}


static void rev_fak(
        const struct kcc_joint *joint,
        const joint_acceleration *qdd,
        struct gc_acc_twist *xdd)
{
    // initialize to 0.0
    memset(xdd->angular_acceleration, 0, sizeof(*xdd->angular_acceleration));
    memset(xdd->linear_acceleration, 0, sizeof(*xdd->linear_acceleration));

    enum joint_axis axis = joint->revolute_joint.axis;
    xdd->angular_acceleration->x = (axis == JOINT_AXIS_X) ? qdd[0] : 0.0;
    xdd->angular_acceleration->y = (axis == JOINT_AXIS_Y) ? qdd[0] : 0.0;
    xdd->angular_acceleration->z = (axis == JOINT_AXIS_Z) ? qdd[0] : 0.0;
}


static void rev_inertial_acceleration(
        const struct kcc_joint *joint,
        const struct gc_twist *xd,
        const joint_velocity *qd,
        struct gc_acc_twist *xdd)
{
    assert(joint);
    assert(xd);
    assert(qd);
    assert(xdd);
    assert(xd->angular_velocity);
    assert(xd->linear_velocity);
    assert(xdd->angular_acceleration);
    assert(xdd->linear_acceleration);

    // Bias acceleration
    //       w_1 x w_2       -> e.g. w_1 x [0, 0, 1]
    // w_1 x v_2 + v_1 x w_2 -> e.g. w_1 x [0, 0, 0] + v_1 x [0, 0, 1] = v_1 x [0, 0, 1]
    if (joint->revolute_joint.axis == JOINT_AXIS_X) {
        double w3m1 = xd->angular_velocity->z * qd[0];
        double w2m1 = xd->angular_velocity->y * qd[0];
        double v3m1 = xd->linear_velocity->z * qd[0];
        double v2m1 = xd->linear_velocity->y * qd[0];

        xdd->angular_acceleration->x =  0.0;
        xdd->angular_acceleration->y =  w3m1;
        xdd->angular_acceleration->z = -w2m1;
        xdd->linear_acceleration->x =  0.0;
        xdd->linear_acceleration->y =  v3m1;
        xdd->linear_acceleration->z = -v2m1;
    } else if (joint->revolute_joint.axis == JOINT_AXIS_Y) {
        double w3m2 = xd->angular_velocity->z * qd[0];
        double w1m2 = xd->angular_velocity->x * qd[0];
        double v3m2 = xd->linear_velocity->z * qd[0];
        double v1m2 = xd->linear_velocity->x * qd[0];

        xdd->angular_acceleration->x = -w3m2;
        xdd->angular_acceleration->y =  0.0;
        xdd->angular_acceleration->z =  w1m2;
        xdd->linear_acceleration->x = -v3m2;
        xdd->linear_acceleration->y =  0.0;
        xdd->linear_acceleration->z =  v1m2;
    } else if (joint->revolute_joint.axis == JOINT_AXIS_Z) {
        double w2m3 = xd->angular_velocity->y * qd[0];
        double w1m3 = xd->angular_velocity->x * qd[0];
        double v2m3 = xd->linear_velocity->y * qd[0];
        double v1m3 = xd->linear_velocity->x * qd[0];

        xdd->angular_acceleration->x =  w2m3;
        xdd->angular_acceleration->y = -w1m3;
        xdd->angular_acceleration->z =  0.0;
        xdd->linear_acceleration->x =  v2m3;
        xdd->linear_acceleration->y = -v1m3;
        xdd->linear_acceleration->z =  0.0;
    } else {
        assert(0);
    }
}


const struct kcc_joint_operators kcc_joint[] = {
    [JOINT_TYPE_REVOLUTE] = {
        .fpk = rev_fpk,
        .fvk = rev_fvk,
        .fak = rev_fak,
        .inertial_acceleration = rev_inertial_acceleration
    }
};