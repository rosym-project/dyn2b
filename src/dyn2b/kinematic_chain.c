#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/linear_algebra.h>
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


void kca_ifk(
        const struct kca_joint *joint,
        const struct ma_wrench *f)
{
    assert(joint);
    assert(f);
    assert(f->frame);
    assert(f->point == f->frame->origin);       // wrench
    assert(joint->target_body == f->body);
    assert(joint->target_frame == f->frame);
}


void kca_ffd(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        struct ma_wrench *f)
{
    assert(joint);
    assert(m);
    assert(f);
    assert(m->frame);
    assert(m->point == m->frame->origin);
    assert(joint->target_body == m->body);
    assert(joint->target_frame == m->frame);

    f->body = joint->reference_body;
    f->frame = joint->target_frame;
    f->point = joint->target_frame->origin;
}


void kca_project_inertia(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        struct ma_abi *r)
{
    assert(joint);
    assert(m);
    assert(r);
    assert(m != r);
    assert(m->frame);
    assert(m->point == m->frame->origin);
    assert(joint->target_body == m->body);
    assert(joint->target_frame == m->frame);

    r->body = joint->reference_body;
    r->point = m->point;
    r->frame = m->frame;
}


void kca_project_wrench(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        const struct ma_wrench *f,
        struct ma_wrench *r)
{
    assert(joint);
    assert(m);
    assert(f);
    assert(r);
    assert(f != r);
    assert(m->frame);
    assert(f->frame);
    assert(m->point == m->frame->origin);
    assert(f->point == f->frame->origin);       // wrench
    assert(joint->target_body == m->body);
    assert(joint->target_frame == m->frame);
    assert(joint->target_body == f->body);
    assert(joint->target_frame == f->frame);

    // Reference: [Featherstone2008]: p. 127, Eq. (7.25)
    r->body = joint->reference_body;
    r->point = m->point;
    r->frame = m->frame;
}


void kca_project_acc_twist(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        const struct ga_acc_twist *xdd,
        struct ga_acc_twist *r)
{
    assert(joint);
    assert(m);
    assert(xdd);
    assert(r);
    assert(xdd != r);
    assert(m->frame);
    assert(xdd->frame);
    assert(m->point == m->frame->origin);
    assert(xdd->point == xdd->frame->origin);       // screw twist
    assert(joint->target_body == m->body);
    assert(joint->target_frame == m->frame);
    assert(joint->target_body == xdd->target_body);
    assert(joint->target_frame == xdd->frame);

    r->target_body = xdd->target_body;
    r->reference_body = xdd->reference_body;
    r->point = xdd->point;
    r->frame = xdd->frame;
}


void kca_cart_force_to_eacc(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        const struct ma_wrench *f1,
        const struct ma_wrench *f2)
{
    assert(joint);
    assert(m);
    assert(f1);
    assert(f2);
    assert(m->frame);
    assert(f1->frame);
    assert(f2->frame);
    assert(m->point == m->frame->origin);
    assert(f1->point == f1->frame->origin);         // wrench
    assert(f2->point == f2->frame->origin);         // wrench
    assert(joint->target_body == m->body);
    assert(joint->target_frame == m->frame);
    assert(joint->target_body == f1->body);
    assert(joint->target_frame == f1->frame);
    assert(joint->target_body == f2->body);
    assert(joint->target_frame == f2->frame);
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


static void rev_ifk(
        const struct kcc_joint *joint,
        const struct mc_wrench *f,
        joint_torque *tau,
        int count)
{
    assert(joint);
    assert(f);
    assert(tau);

    int k = joint->revolute_joint.axis;

    for (int i = 0; i < count; i++) {
        tau[i] = f->torque[i].data[k];
    }
}


static void rev_ffd(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        const joint_torque *tau,
        struct mc_wrench *f,
        int count)
{
    assert(joint);
    assert(m);
    assert(tau);
    assert(f);

    int k = joint->revolute_joint.axis;
    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];

    for (int i = 0; i < count; i++) {
        double qdd = tau[i] / d;
        for (int j = 0; j < 3; j++) {
            f->torque[i].data[j] = m->second_moment_of_mass.row[j].data[k] * qdd;
            f->force[i].data[j] = m->first_moment_of_mass.row[k].data[j] * qdd;    // consider transpose, thus [k, j]
        }
    }
}


static void rev_project_inertia(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        struct mc_abi *r)
{
    assert(joint);
    assert(m);
    assert(r);
    assert(m != r);

    int k = joint->revolute_joint.axis;

    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];
    assert(d != 0.0);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            // The (__ik/d) element represents an entry in the projection matrix P
            // The __ij and __kj elements represent the entries of the inertia matrix M
            //     __ij is associated with the diagonal one-element of the projection matrix (the "identity" matrix part)
            //     -__kj is associated with the non-zero column of the projection matrix

            // 0th moment of mass matrix
            double m0ij = m->zeroth_moment_of_mass.row[i].data[j];
            double m0ik = m->first_moment_of_mass.row[k].data[i];    // consider transpose, thus [k, i]
            double m0kj = m->first_moment_of_mass.row[k].data[j];

            double pm0 = m0ij - (m0ik * m0kj) / d;

            r->zeroth_moment_of_mass.row[i].data[j] = pm0;


            // 1st moment of mass matrix
            double m1ij = m->first_moment_of_mass.row[i].data[j];
            double m1ik = m->second_moment_of_mass.row[i].data[k];
            double m1kj = m->first_moment_of_mass.row[k].data[j];
            double pm1 = m1ij - (m1ik * m1kj) / d;

            r->first_moment_of_mass.row[i].data[j] = pm1;


            // 2nd moment of mass matrix
            double mij = m->second_moment_of_mass.row[i].data[j];
            double mik = m->second_moment_of_mass.row[i].data[k];
            double mkj = m->second_moment_of_mass.row[k].data[j];
            double pm = mij - (mik * mkj) / d;

            r->second_moment_of_mass.row[i].data[j] = pm;
        }
    }
}


static void rev_project_wrench(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        const struct mc_wrench *f,
        struct mc_wrench *r,
        int count)
{
    assert(joint);
    assert(m);
    assert(f);
    assert(r);
    assert(f != r);

    int k = joint->revolute_joint.axis;

    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];
    assert(d != 0.0);

    for (int j = 0; j < count; j++) {
        double qdd = f->torque[j].data[k] / d;

        for (int i = 0; i < 3; i++) {
            double m2 = m->second_moment_of_mass.row[i].data[k];
            r->torque[j].data[i] = f->torque[j].data[i] - (m2 * qdd);

            double m1 = m->first_moment_of_mass.row[k].data[i];               // consider transpose, thus [k, i]
            r->force[j].data[i] = f->force[j].data[i]- (m1 * qdd);
        }
    }
}


static void rev_project_acc_twist(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        const struct gc_acc_twist *xdd,
        struct gc_acc_twist *r)
{
    assert(joint);
    assert(m);
    assert(xdd);
    assert(r);
    assert(xdd != r);

    int k = joint->revolute_joint.axis;

    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];
    assert(d != 0.0);

    double tau = 0.0;

    for (int i = 0; i < 3; i++) {
        double m2 = m->second_moment_of_mass.row[i].data[k];
        double wd = xdd->angular_acceleration[0].data[i];
        tau += m2 * wd;
        r->angular_acceleration[0].data[i] = wd;

        double m1 = m->first_moment_of_mass.row[i].data[k];
        double vd = xdd->linear_acceleration[0].data[i];
        tau += m1 * vd;
        r->linear_acceleration[0].data[i] = vd;
    }

    r->angular_acceleration[0].data[k] -= tau / d;
}


static void rev_cart_force_to_eacc(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        const struct mc_wrench *f1,
        const struct mc_wrench *f2,
        mc_eacc *e,
        int count_f1,
        int count_f2)
{
    assert(joint);
    assert(m);
    assert(f1);
    assert(f2);
    assert(e);

    int k = joint->revolute_joint.axis;

    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];
    assert(d != 0.0);

    la_dger_os(count_f1, count_f2,
            1.0 / d,
            &f1->torque[0].data[k], 3,
            &f2->torque[0].data[k], 3,
            e, count_f2);
}


static void rev_decomp_e_cstr(
        const struct kcc_joint *joint,
        const struct mc_abi *m,
        const struct mc_wrench *f,
        double *d_cstr,
        double *r,
        int count)
{
    assert(joint);
    assert(m);
    assert(f);
    assert(r);

    int k = joint->revolute_joint.axis;

    double d = m->second_moment_of_mass.row[k].data[k]
            + joint->revolute_joint.inertia[0];
    assert(d != 0.0);

    la_dsytrfr_lo(count,
            1.0 / d, &f->torque[0].data[k], 3,
            d_cstr, count,
            r, count);
}


const struct kcc_joint_operators kcc_joint[] = {
    [JOINT_TYPE_REVOLUTE] = {
        .fpk = rev_fpk,
        .fvk = rev_fvk,
        .fak = rev_fak,
        .inertial_acceleration = rev_inertial_acceleration,
        .ifk = rev_ifk,
        .ffd = rev_ffd,
        .project_inertia = rev_project_inertia,
        .project_wrench = rev_project_wrench,
        .project_acc_twist = rev_project_acc_twist,
        .cart_force_to_eacc = rev_cart_force_to_eacc,
        .decomp_e_cstr = rev_decomp_e_cstr
    }
};
