#include <dyn2b/functions/mechanics.h>
#include <dyn2b/functions/linear_algebra.h>

#include <math.h>
#include <stdio.h>
#include <assert.h>


void mc_eacc_add(
        const mc_eacc *e1,
        const mc_eacc *e2,
        mc_eacc *e,
        int m,
        int n)
{
    assert(e1);
    assert(e2);
    assert(e);

    la_dgeadd_os(m, n,
            e1, n,
            e2, n,
            e, n);
}


void mc_eacc_sub(
        const mc_eacc *e1,
        const mc_eacc *e2,
        mc_eacc *e,
        int m,
        int n)
{
    assert(e1);
    assert(e2);
    assert(e);

    la_dgeadd_oe(m, n, -1.0,
            e2, n,
            e1, n,
            e, n);
}


void mc_eacc_balance(
        const double *decomp,
        const mc_eacc *e,
        double *scale,
        int count)
{
    double z[count];
    la_trsv_lnd(count,
            decomp, count,
            e, 1,
            z, 1);

    for (int i = 0; i < count; i++) {
        const double EPSILON = 1e-15;
        double d = decomp[i * count + i];

        // Truncation
        if (fabs(d) < EPSILON) z[i] = 0.0;
        else z[i] /= d;
    }

    la_trsv_ltd(count,
            decomp, count,
            z, 1,
            scale, 1);
}


void mc_momentum_derive(
        const struct gc_twist *xd,
        const struct mc_momentum *p,
        struct mc_wrench *r)
{
    // w x n
    struct vector3 wxn;
    la_dcross_o(
            (double *)xd->angular_velocity, 1,
            (double *)p->angular_momentum, 1,
            (double *)&wxn, 1);

    // v x f
    struct vector3 vxf;
    la_dcross_o(
            (double *)xd->linear_velocity, 1,
            (double *)p->linear_momentum, 1,
            (double *)&vxf, 1);

    // n' = w x n + v x f
    la_daxpy_oe(3,
            1.0, (double *)&wxn, 1,
            (double *)&vxf, 1,
            (double *)r->torque, 1);

    // f' = w x f
    la_dcross_o(
            (double *)xd->angular_velocity, 1,
            (double *)p->linear_momentum, 1,
            (double *)r->force, 1);
}


void ma_momentum_derive(
        const struct ga_twist *xd,
        const struct ma_momentum *p,
        struct ma_wrench *r)
{
    assert(xd);
    assert(p);
    assert(r);
    assert(xd->frame);
    assert(p->frame);
    assert(xd->frame->origin == xd->point);     // screw twist
    assert(p->frame->origin == p->point);       // wrench
    assert(xd->target_body == p->body);
    assert(xd->point == p->point);
    assert(xd->frame == p->frame);

    r->body = xd->target_body;
    r->point = xd->point;
    r->frame = xd->frame;
}


void mc_wrench_tf_tgt_to_ref(
        const struct gc_pose *x,
        const struct mc_wrench *f,
        struct mc_wrench *r,
        int count)
{
    assert(x);
    assert(f);
    assert(r);
    assert(f != r);
    assert(x->rotation && x->translation);
    assert(f->torque && f->force);
    assert(r->torque && r->force);
    assert(count >= 0);

    // f' = E^T f
    la_dgemm_nnos(count, 3, 3,
            (double *)f->force, 3,
            (double *)x->rotation, 3,
            (double *)r->force, 3);

    // rx E^T f = r x f'
    struct vector3 rxetf[count];
    for (int i = 0; i < count; i++) {
        la_dcross_o(
                (double *)x->translation, 1,
                (double *)&r->force[i], 1,
                (double *)&rxetf[i], 1);
    }

    // n' = E^T n + rx E^T f
    la_dgemm_nnoe(count, 3, 3,
            1.0, (double *)f->torque, 3, (double *)x->rotation, 3,
            1.0, (double *)&rxetf[0], 3,
            (double *)r->torque, 3);
}


void ma_wrench_tf_tgt_to_ref(
        const struct ga_pose *x,
        const struct ma_wrench *f,
        struct ma_wrench *r)
{
    assert(x);
    assert(f);
    assert(r);
    assert(f != r);
    assert(f->point == f->frame->origin);
    assert(f->frame == x->target_frame);

    r->body = f->body;
    r->frame = x->reference_frame;
    r->point = r->frame->origin;
}


void mc_wrench_invert(
        const struct mc_wrench *f,
        struct mc_wrench *r,
        int count)
{
    assert(f);
    assert(r);

    la_dscal_o(3 * count, -1.0, (double *)f->torque, 1, (double *)r->torque, 1);
    la_dscal_o(3 * count, -1.0, (double *)f->force, 1, (double *)r->force, 1);
}


void ma_wrench_invert(
        const struct ma_wrench *f,
        struct ma_wrench *r)
{
    assert(f);
    assert(r);
    assert(f->frame);
    assert(f->frame->origin == f->point);

    r->body = f->body;
    r->point = f->point;
    r->frame = f->frame;
}


void mc_wrench_add(
        const struct mc_wrench *f1,
        const struct mc_wrench *f2,
        struct mc_wrench *r,
        int count)
{
    assert(f1);
    assert(f2);
    assert(r);

    la_daxpy_oe(3 * count,
            1.0, (double *)f1->torque, 1, (double *)f2->torque, 1,
            (double *)r->torque, 1);
    la_daxpy_oe(3 * count,
            1.0, (double *)f1->force, 1, (double *)f2->force, 1,
            (double *)r->force, 1);
}


void ma_wrench_add(
        const struct ma_wrench *f1,
        const struct ma_wrench *f2,
        struct ma_wrench *r)
{
    assert(f1);
    assert(f2);
    assert(r);
    assert(f1->frame);
    assert(f2->frame);
    assert(f1->frame->origin == f1->point);
    assert(f2->frame->origin == f2->point);
    assert(f1->body == f2->body);
    assert(f1->point == f2->point);
    assert(f1->frame == f2->frame);

    r->body = f1->body;
    r->point = f1->point;
    r->frame = f1->frame;
}


void mc_wrench_sub(
        const struct mc_wrench *f1,
        const struct mc_wrench *f2,
        struct mc_wrench *r,
        int count)
{
    assert(f1);
    assert(f2);
    assert(r);

    la_daxpy_oe(3 * count,
            -1.0, (double *)f2->torque, 1, (double *)f1->torque, 1,
            (double *)r->torque, 1);
    la_daxpy_oe(3 * count,
            -1.0, (double *)f2->force, 1, (double *)f1->force, 1,
            (double *)r->force, 1);
}


void ma_wrench_sub(
        const struct ma_wrench *f1,
        const struct ma_wrench *f2,
        struct ma_wrench *r)
{
    assert(f1);
    assert(f2);
    assert(r);
    assert(f1->frame);
    assert(f2->frame);
    assert(f1->frame->origin == f1->point);
    assert(f2->frame->origin == f2->point);
    assert(f1->body == f2->body);
    assert(f1->point == f2->point);
    assert(f1->frame == f2->frame);

    r->body = f1->body;
    r->point = f1->point;
    r->frame = f1->frame;
}


void mc_wrench_log(
        const struct mc_wrench *f,
        int count)
{
    assert(f);
    assert(f->torque);
    assert(f->force);

    printf("WrenchCoord(\n");
    for (int i = 0; i < count; i++) {
        printf("  torque=[%5.2f, %5.2f, %5.2f], ",
                f->torque[i].x,
                f->torque[i].y,
                f->torque[i].z);
        printf("force=[%5.2f, %5.2f, %5.2f]",
                f->force[i].x,
                f->force[i].y,
                f->force[i].z);
        if (i != count - 1) printf("\n");
    }
    printf(")\n");
}


void ma_wrench_log(
        const struct ma_wrench *f)
{
    assert(f);
    assert(f->point);
    assert(f->body);
    assert(f->frame);

    printf("WrenchADT(target=%s|%s, frame={%s})\n",
        f->point->name,
        f->body->name,
        f->frame->name);
}


void mc_rbi_map_twist_to_momentum(
        const struct mc_rbi *m,
        const struct gc_twist *xd,
        struct mc_momentum *r)
{
    assert(m);
    assert(xd);
    assert(r);

    // h x v
    struct vector3 hxv;
    la_dcross_o(
            (double *)&m->first_moment_of_mass, 1,
            (double *)xd->linear_velocity, 1,
            (double *)&hxv, 1);

    // n = I w + h x v
    la_dgemv_noe(3, 3,
            1.0, (double *)&m->second_moment_of_mass, 3, (double *)xd->angular_velocity, 1,
            1.0, (double *)&hxv, 1,
            (double *)r->angular_momentum, 1);

    // h x w
    struct vector3 hxw;
    la_dcross_o(
            (double *)&m->first_moment_of_mass, 1,
            (double *)xd->angular_velocity, 1,
            (double *)&hxw, 1);

    // m v
    la_dscal_o(3,
            m->zeroth_moment_of_mass,
            (double *)xd->linear_velocity, 1,
            (double *)r->linear_momentum, 1);

    // f = m v - h x w
    la_daxpy_ie(3,
            -1.0, (double *)&hxw, 1,
            (double *)r->linear_momentum, 1);
}


void ma_rbi_map_twist_to_momentum(
        const struct ma_rbi *m,
        const struct ga_twist *xd,
        struct ma_momentum *r)
{
    assert(m);
    assert(xd);
    assert(r);
    assert(xd->frame);
    assert(xd->frame->origin == xd->point);       // screw twist
    assert(m->body == xd->target_body);
    assert(m->point == xd->point);
    assert(m->frame == xd->frame);

    r->body = xd->target_body;
    r->point = xd->point;
    r->frame = xd->frame;
}


void mc_rbi_to_abi(
        const struct mc_rbi *rbi,
        struct mc_abi *r)
{
    assert(rbi);
    assert(r);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            r->zeroth_moment_of_mass.row[i].data[j] = 0.0;
            r->first_moment_of_mass.row[i].data[j] = 0.0;
            r->second_moment_of_mass.row[i].data[j]
                    = rbi->second_moment_of_mass.row[i].data[j];
        }
        r->zeroth_moment_of_mass.row[i].data[i] = rbi->zeroth_moment_of_mass;
    }
    la_dcrossop(
            (double *)&rbi->first_moment_of_mass, 1,
            (double *)&r->first_moment_of_mass, 3);
}


void ma_rbi_to_abi(
        const struct ma_rbi *rbi,
        struct ma_abi *r)
{
    assert(rbi);
    assert(r);

    r->body = rbi->body;
    r->point = rbi->point;
    r->frame = rbi->frame;
}


void mc_rbi_log(
        const struct mc_rbi *m)
{
    assert(m);

    printf("RigidBodyInertiaCoord(\n");
    printf("  I=                     h=       m=\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%5.2f, %5.2f, %5.2f]",
                m->second_moment_of_mass.row[i].x,
                m->second_moment_of_mass.row[i].y,
                m->second_moment_of_mass.row[i].z);
        printf("  [%5.2f]", m->first_moment_of_mass.data[i]);
        if (i == 0) printf("  [%5.2f]", m->zeroth_moment_of_mass);
        if (i != 2) printf(",\n");
    }
    printf(")\n");
}


void ma_rbi_log(
        const struct ma_rbi *m)
{
    assert(m);
    assert(m->body);
    assert(m->point);
    assert(m->frame);

    printf("RigidBodyInertiaADT(point=%s|%s, frame={%s})\n",
        m->point->name,
        m->body->name,
        m->frame->name);
}


void mc_abi_tf_tgt_to_ref(
        const struct gc_pose *x,
        const struct mc_abi *m,
        struct mc_abi *r)
{
    assert(x);
    assert(m);
    assert(r);
    assert(m != r);

    // M' = E^T M E
    struct matrix3x3 me;
    la_dgemm_nnos(3, 3, 3,
            (double *)&m->zeroth_moment_of_mass, 3,
            (double *)x->rotation, 3,
            (double *)&me, 3);
    la_dgemm_tnos(3, 3, 3,
            (double *)x->rotation, 3,
            (double *)&me, 3,
            (double *)&r->zeroth_moment_of_mass, 3);

    // H' = E^T H E + rxM'
    struct matrix3x3 he;
    struct matrix3x3 ethe;
    struct matrix3x3 rx;

    la_dgemm_nnos(3, 3, 3,
            (double *)&m->first_moment_of_mass, 3,
            (double *)x->rotation, 3,
            (double *)&he, 3);
    la_dgemm_tnos(3, 3, 3,
            (double *)x->rotation, 3,
            (double *)&he, 3,
            (double *)&ethe, 3);
    la_dcrossop(
            (double *)x->translation, 1,
            (double *)&rx, 3);
    la_dgemm_nnoe(3, 3, 3,
            1.0, (double *)&rx, 3, (double *)&r->zeroth_moment_of_mass, 3,
            1.0, (double *)&ethe, 3,
            (double *)&r->first_moment_of_mass, 3);

    // I' = E^T I E + rx(E^T H E)^T - H'rx
    struct matrix3x3 ie;
    struct matrix3x3 etie;
    struct matrix3x3 etie_rxethet;

    // E^T I E
    la_dgemm_nnos(3, 3, 3,
            (double *)&m->second_moment_of_mass, 3,
            (double *)x->rotation, 3,
            (double *)&ie, 3);
    la_dgemm_tnos(3, 3, 3,
            (double *)x->rotation, 3,
            (double *)&ie, 3,
            (double *)&etie, 3);

    // + rx(E^T H E)^T
    la_dgemm_ntoe(3, 3, 3,
            1.0, (double *)&rx, 3, (double *)&ethe, 3,
            1.0, (double *)&etie, 3,
            (double *)&etie_rxethet, 3);

    // - H'rx
    la_dgemm_nnoe(3, 3, 3,
            -1.0, (double *)&r->first_moment_of_mass, 3, (double *)&rx, 3,
            1.0, (double *)&etie_rxethet, 3,
            (double *)&r->second_moment_of_mass, 3);
}


void ma_abi_tf_tgt_to_ref(
        const struct ga_pose *x,
        const struct ma_abi *m,
        struct ma_abi *r)
{
    assert(x);
    assert(m);
    assert(r);
    assert(m != r);
    assert(m->point == m->frame->origin);
    assert(m->frame == x->target_frame);

    r->body = m->body;
    r->frame = x->reference_frame;
    r->point = r->frame->origin;
}


void mc_abi_add(
        const struct mc_abi *m1,
        const struct mc_abi *m2,
        struct mc_abi *r)
{
    assert(m1);
    assert(m2);
    assert(r);

    la_dgeadd_os(3, 3,
            (double *)&m1->zeroth_moment_of_mass, 3,
            (double *)&m2->zeroth_moment_of_mass, 3,
            (double *)&r->zeroth_moment_of_mass, 3);
    la_dgeadd_os(3, 3,
            (double *)&m1->first_moment_of_mass, 3,
            (double *)&m2->first_moment_of_mass, 3,
            (double *)&r->first_moment_of_mass, 3);
    la_dgeadd_os(3, 3,
            (double *)&m1->second_moment_of_mass, 3,
            (double *)&m2->second_moment_of_mass, 3,
            (double *)&r->second_moment_of_mass, 3);
}


void ma_abi_add(
        const struct ma_abi *m1,
        const struct ma_abi *m2,
        struct ma_abi *r)
{
    assert(m1);
    assert(m2);
    assert(r);
    assert(m1->body == m2->body);
    assert(m1->point == m2->point);
    assert(m1->frame == m2->frame);

    r->body = m1->body;
    r->point = m1->point;
    r->frame = m1->frame;
}


void mc_abi_map_acc_twist_to_wrench(
        const struct mc_abi *m,
        const struct gc_acc_twist *xdd,
        struct mc_wrench *f)
{
    assert(m);
    assert(xdd);
    assert(f);

    // H v
    struct vector3 hv;
    la_dgemv_nos(3, 3,
            (double *)&m->first_moment_of_mass, 3,
            (double *)xdd->linear_acceleration, 1,
            (double *)&hv, 1);

    // n = I w + H v
    la_dgemv_noe(3, 3,
            1.0, (double *)&m->second_moment_of_mass, 3, (double *)xdd->angular_acceleration, 1,
            1.0, (double *)&hv, 1,
            (double *)f->torque, 1);

    // H^T w
    struct vector3 htw;
    la_dgemv_tos(3, 3,
            (double *)&m->first_moment_of_mass, 3,
            (double *)xdd->angular_acceleration, 1,
            (double *)&htw, 1);

    // f = M v + H^T w
    la_dgemv_noe(3, 3,
            1.0, (double *)&m->zeroth_moment_of_mass, 3, (double *)xdd->linear_acceleration, 1,
            1.0, (double *)&htw, 1,
            (double *)f->force, 1);
}


void ma_abi_map_acc_twist_to_wrench(
        const struct ma_abi *m,
        const struct ga_acc_twist *xdd,
        struct ma_wrench *f)
{
    assert(m);
    assert(xdd);
    assert(f);
    assert(xdd->frame);
    assert(xdd->point == xdd->frame->origin);       // screw twist
    assert(xdd->point == m->point);
    assert(xdd->target_body == m->body);
    assert(xdd->frame == m->frame);

    f->body = xdd->target_body;
    f->point = xdd->point;
    f->frame = xdd->frame;
}


void mc_abi_log(
        const struct mc_abi *m)
{
    assert(m);

    printf("ArticulatedBodyInertiaCoord(\n");
    printf("  I=                     H=                     M=\n");
    for (int i = 0; i < 3; i++) {
        printf("  [%5.2f, %5.2f, %5.2f]",
                m->second_moment_of_mass.row[i].x,
                m->second_moment_of_mass.row[i].y,
                m->second_moment_of_mass.row[i].z);
        printf("  [%5.2f, %5.2f, %5.2f]",
                m->first_moment_of_mass.row[i].x,
                m->first_moment_of_mass.row[i].y,
                m->first_moment_of_mass.row[i].z);
        printf("  [%5.2f, %5.2f, %5.2f]",
                m->zeroth_moment_of_mass.row[i].x,
                m->zeroth_moment_of_mass.row[i].y,
                m->zeroth_moment_of_mass.row[i].z);
        if (i != 2) printf(",\n");
    }
    printf(")\n");
}


void ma_abi_log(
        const struct ma_abi *m)
{
    assert(m);
    assert(m->body);
    assert(m->point);
    assert(m->frame);

    printf("ArticulatedBodyInertiaADT(point=%s|%s, frame={%s})\n",
        m->point->name,
        m->body->name,
        m->frame->name);
}
