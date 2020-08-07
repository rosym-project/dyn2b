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
