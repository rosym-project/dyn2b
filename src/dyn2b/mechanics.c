#include <dyn2b/functions/mechanics.h>
#include <dyn2b/functions/linear_algebra.h>

#include <stdio.h>
#include <assert.h>


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