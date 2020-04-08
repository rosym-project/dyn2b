#include <dyn2b/example/solver_state.h>

#include <stdlib.h>


/**
 * Setup the solver state for a simple task and a serial kinematic chain
 * (coordinates).
 */
void setup_simple_state_c(
        const struct kcc_kinematic_chain *kc,
        struct solver_state_c *s)
{
    const int NR_SEGMENTS = kc->number_of_segments;
    const int NR_SEGMENTS_WITH_BASE = NR_SEGMENTS + 1;
    const int EE = NR_SEGMENTS;

    // FPK
    s->nbody = NR_SEGMENTS;
    s->nq = NR_SEGMENTS;

    // FPK
    s->x_jnt = calloc(NR_SEGMENTS, sizeof(struct gc_pose));
    s->x_rel = calloc(NR_SEGMENTS, sizeof(struct gc_pose));
    s->x_tot = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct gc_pose));

    // FPK
    s->q     = calloc(s->nq, sizeof(double));

    for (int i = 0; i < NR_SEGMENTS; i++) {
        // FPK
        s->x_jnt[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_jnt[i].translation = calloc(1, sizeof(struct vector3));
        s->x_rel[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_rel[i].translation = calloc(1, sizeof(struct vector3));
    }

    for (int i = 0; i < NR_SEGMENTS_WITH_BASE; i++) {
        // FPK
        s->x_tot[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_tot[i].translation = calloc(1, sizeof(struct vector3));
    }

    // FPK
    s->x_tot[0].rotation->row_x.x = 1.0;
    s->x_tot[0].rotation->row_y.y = 1.0;
    s->x_tot[0].rotation->row_z.z = 1.0;
}


/**
 * Setup the solver state for a simple task and a serial kinematic chain (ADT).
 */
void setup_simple_state_a(
        const struct kca_kinematic_chain *kc,
        struct solver_state_a *s)
{
    const int NR_SEGMENTS = kc->number_of_segments;
    const int NR_SEGMENTS_WITH_BASE = NR_SEGMENTS + 1;
    const int EE = NR_SEGMENTS;

    // FPK
    s->nbody = NR_SEGMENTS;

    // FPK
    s->x_jnt = calloc(NR_SEGMENTS, sizeof(struct ga_pose));
    s->x_rel = calloc(NR_SEGMENTS, sizeof(struct ga_pose));
    s->x_tot = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ga_pose));

    struct body *body_world = kc->segment[0].joint_attachment.reference_body;
    struct frame *frame_world = kc->segment[0].joint_attachment.reference_frame;
    struct point *point_world_origin = frame_world->origin;
    struct body *body_ee = kc->segment[NR_SEGMENTS - 1].joint.target_body;
    struct frame *frame_ee = kc->segment[NR_SEGMENTS - 1].joint.target_frame;
    struct point *point_ee_origin = frame_ee->origin;

    // FPK
    s->x_tot[0].target_body = body_world;
    s->x_tot[0].reference_body = body_world;
    s->x_tot[0].target_frame = frame_world;
    s->x_tot[0].reference_frame = frame_world;
}