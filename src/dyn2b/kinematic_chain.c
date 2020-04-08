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


const struct kcc_joint_operators kcc_joint[] = {
    [JOINT_TYPE_REVOLUTE] = {
        .fpk = rev_fpk,
        .fvk = rev_fvk
    }
};