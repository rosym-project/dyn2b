#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/linear_algebra.h>

#include <stdio.h>
#include <assert.h>


void gc_pose_compose(
        const struct gc_pose *x1,
        const struct gc_pose *x2,
        struct gc_pose *r)
{
    assert(x1);
    assert(x2);
    assert(r);
    assert(x1 != r);
    assert(x2 != r);
    assert(x1->rotation && x1->translation);
    assert(x2->rotation && x2->translation);
    assert(r->rotation && r->translation);

    // E' = E_1 E_2
    la_dgemm_nnos(3, 3, 3,
            (double *)x1->rotation, 3,
            (double *)x2->rotation, 3,
            (double *)r->rotation, 3);

    // r' = r_2 + E_2^T r_1
    la_dgemv_toe(3, 3,
            1.0, (double *)x2->rotation, 3, (double *)x1->translation, 1,
            1.0, (double *)x2->translation, 1,
            (double *)r->translation, 1);
}


void ga_pose_compose(
        const struct ga_pose *x1,
        const struct ga_pose *x2,
        struct ga_pose *r)
{
    assert(x1);
    assert(x2);
    assert(r);
    assert(x1 != r);
    assert(x2 != r);
    assert(x1->reference_frame == x2->target_frame);

    r->target_body = x1->target_body;
    r->target_frame = x1->target_frame;
    r->reference_body = x2->reference_body;
    r->reference_frame = x2->reference_frame;
}


void gc_pose_log(
        const struct gc_pose *x)
{
    assert(x);
    assert(x->rotation && x->translation);

    printf("PoseCoord(rotation=[");
    for (int i = 0; i < 3; i++) {
        if (i != 0) printf("                    ");
        printf("[%5.2f, %5.2f, %5.2f]",
                x->rotation->row[i].x,
                x->rotation->row[i].y,
                x->rotation->row[i].z);
        if (i != 2) printf(",\n");
    }
    printf("],\n       translation=[[%5.3f, %5.3f, %5.3f]])\n",
            x->translation->x,
            x->translation->y,
            x->translation->z);
}


void ga_pose_log(
        const struct ga_pose *x)
{
    assert(x);
    assert(x->target_frame);
    assert(x->target_body);
    assert(x->reference_frame);
    assert(x->reference_body);

    printf("PoseADT(target={%s}|%s, reference={%s}|%s)\n",
            x->target_frame->name,
            x->target_body->name,
            x->reference_frame->name,
            x->reference_body->name);
}


void gc_twist_tf_ref_to_tgt(
        const struct gc_pose *x,
        const struct gc_twist *xd,
        struct gc_twist *r)
{
    assert(x);
    assert(xd);
    assert(r);
    assert(xd != r);

    // r x w
    struct vector3 rxw;
    la_dcross_o((double *)x->translation, 1,
            (double *)xd->angular_velocity, 1,
            (double *)&rxw, 1);

    // v - r x w
    struct vector3 v_rxw;
    la_daxpy_oe(3,
            -1.0, (double *)&rxw, 1,
            (double *)xd->linear_velocity, 1,
            (double *)&v_rxw, 1);

    // v' = E(v - r x w)
    la_dgemv_nos(3, 3,
            (double *)x->rotation, 3,
            (double *)xd->angular_velocity, 1,
            (double *)r->angular_velocity, 1);

    // w' = E w
    la_dgemv_nos(3, 3,
            (double *)x->rotation, 3,
            (double *)&v_rxw, 1,
            (double *)r->linear_velocity, 1);
}


void ga_twist_tf_ref_to_tgt(
        const struct ga_pose *x,
        const struct ga_twist *xd,
        struct ga_twist *r)
{
    assert(x);
    assert(xd);
    assert(r);
    assert(xd != r);
    assert(xd->point == xd->frame->origin);       // screw twist
    assert(xd->frame == x->reference_frame);

    r->target_body = xd->target_body;
    r->reference_body = xd->reference_body;
    r->frame = x->target_frame;
    r->point = r->frame->origin;
}


void gc_twist_accumulate(
        const struct gc_twist *xd1,
        const struct gc_twist *xd2,
        struct gc_twist *r)
{
    assert(xd1);
    assert(xd2);
    assert(r);

    la_daxpy_oe(3,
            1.0, (double *)xd1->angular_velocity, 1,
            (double *)xd2->angular_velocity, 1,
            (double *)r->angular_velocity, 1);

    la_daxpy_oe(3,
            1.0, (double *)xd1->linear_velocity, 1,
            (double *)xd2->linear_velocity, 1,
            (double *)r->linear_velocity, 1);
}


void ga_twist_accumulate(
        const struct ga_twist *xd1,
        const struct ga_twist *xd2,
        struct ga_twist *r)
{
    assert(xd1);
    assert(xd2);
    assert(r);
    assert(xd1->frame);
    assert(xd2->frame);
    assert(xd1->frame->origin);
    assert(xd2->frame->origin);
    assert(xd1->point == xd1->frame->origin);       // screw twist
    assert(xd2->point == xd2->frame->origin);       // screw twist
    assert(xd1->target_body == xd2->reference_body);
    assert(xd1->point == xd2->point);
    assert(xd1->frame == xd2->frame);

    r->target_body = xd2->target_body;
    r->reference_body = xd1->reference_body;
    r->point = xd1->point;
    r->frame = xd1->frame;
}


void gc_twist_log(
        const struct gc_twist *xd)
{
    assert(xd);
    assert(xd->angular_velocity && xd->linear_velocity);

    printf("TwistCoord(angular=[%5.2f, %5.2f, %5.2f],\n",
            xd->angular_velocity->x,
            xd->angular_velocity->y,
            xd->angular_velocity->z);
    printf("           linear =[%5.2f, %5.2f, %5.2f])\n",
            xd->linear_velocity->x,
            xd->linear_velocity->y,
            xd->linear_velocity->z);
}


void ga_twist_log(
        const struct ga_twist *xd)
{
    assert(xd);
    assert(xd->point);
    assert(xd->target_body);
    assert(xd->reference_body);
    assert(xd->frame);

    printf("TwistADT(target=%s|%s, reference=%s, frame={%s})\n",
            xd->point->name,
            xd->target_body->name,
            xd->reference_body->name,
            xd->frame->name);
}