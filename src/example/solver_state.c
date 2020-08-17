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
    const int NR_CSTR = 6;

    // FPK
    s->nbody = NR_SEGMENTS;
    s->nq = NR_SEGMENTS;
    // FVK
    s->nd = NR_SEGMENTS;

    // FPK
    s->x_jnt = calloc(NR_SEGMENTS, sizeof(struct gc_pose));
    s->x_rel = calloc(NR_SEGMENTS, sizeof(struct gc_pose));
    s->x_tot = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct gc_pose));
    // FVK
    s->xd_jnt = calloc(NR_SEGMENTS, sizeof(struct gc_twist));
    s->xd_tf  = calloc(NR_SEGMENTS, sizeof(struct gc_twist));
    s->xd     = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct gc_twist));
    // FAK
    s->xdd_jnt   = calloc(NR_SEGMENTS, sizeof(struct gc_acc_twist));
    s->xdd_net   = calloc(NR_SEGMENTS, sizeof(struct gc_acc_twist));
    s->xdd_bias  = calloc(NR_SEGMENTS, sizeof(struct gc_acc_twist));
    s->xdd_tf    = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct gc_acc_twist));
    s->xdd_nact  = calloc(NR_SEGMENTS, sizeof(struct gc_acc_twist));
    s->xdd       = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct gc_acc_twist));
    // Inertia
    s->m_art = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct mc_abi));
    s->m_app = calloc(NR_SEGMENTS, sizeof(struct mc_abi));
    s->m_tf  = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct mc_abi));
    // Inertial force
    s->p            = calloc(NR_SEGMENTS, sizeof(struct mc_momentum));
    s->f_bias_art   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct mc_wrench));
    s->f_bias_eom   = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_bias_app   = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_bias_tf    = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_bias_nact  = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->tau_bias_art = calloc(s->nd, sizeof(joint_torque));
    // Feed-forward torque
    s->f_ff_art   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct mc_wrench));
    s->f_ff_app   = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_ff_jnt   = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->tau_ff_art = calloc(s->nd, sizeof(joint_torque));
    // External force
    s->tau_ff      = calloc(s->nd, sizeof(joint_torque));
    s->f_ext       = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_ext_art   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct mc_wrench));
    s->f_ext_app   = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->f_ext_tf    = calloc(NR_SEGMENTS, sizeof(struct mc_wrench));
    s->tau_ext_art = calloc(s->nd, sizeof(joint_torque));
    // Constraint force
    s->f_cstr_art   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_wrench));
    s->f_cstr_app   = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->d_cstr_art   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(double *));
    s->nu_cstr      = calloc(NR_CSTR, sizeof(double));
    s->e_cstr       = calloc(NR_CSTR, sizeof(mc_eacc));
    s->tau_cstr_art = calloc(s->nd, sizeof(joint_torque *));
    s->tau_cstr     = calloc(s->nd, sizeof(joint_torque));

    // FPK
    s->q     = calloc(s->nq, sizeof(double));
    // FVK
    s->qd    = calloc(s->nd, sizeof(double));
    // FAK
    s->qdd   = calloc(s->nd, sizeof(double));
    // Dynamics
    s->tau_ctrl = calloc(s->nd, sizeof(joint_torque));

    for (int i = 0; i < NR_SEGMENTS; i++) {
        // FPK
        s->x_jnt[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_jnt[i].translation = calloc(1, sizeof(struct vector3));
        s->x_rel[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_rel[i].translation = calloc(1, sizeof(struct vector3));
        // FVK
        s->xd_jnt[i].angular_velocity = calloc(1, sizeof(struct vector3));
        s->xd_jnt[i].linear_velocity  = calloc(1, sizeof(struct vector3));
        s->xd_tf[i].angular_velocity  = calloc(1, sizeof(struct vector3));
        s->xd_tf[i].linear_velocity   = calloc(1, sizeof(struct vector3));
        // FAK
        s->xdd_jnt[i].angular_acceleration  = calloc(1, sizeof(struct vector3));
        s->xdd_jnt[i].linear_acceleration   = calloc(1, sizeof(struct vector3));
        s->xdd_bias[i].angular_acceleration = calloc(1, sizeof(struct vector3));
        s->xdd_bias[i].linear_acceleration  = calloc(1, sizeof(struct vector3));
        s->xdd_net[i].angular_acceleration  = calloc(1, sizeof(struct vector3));
        s->xdd_net[i].linear_acceleration   = calloc(1, sizeof(struct vector3));
        s->xdd_nact[i].angular_acceleration = calloc(1, sizeof(struct vector3));
        s->xdd_nact[i].linear_acceleration  = calloc(1, sizeof(struct vector3));
        // Inertial force
        s->p[i].angular_momentum = calloc(1, sizeof(struct vector3));
        s->p[i].linear_momentum  = calloc(1, sizeof(struct vector3));
        s->f_bias_eom[i].torque  = calloc(1, sizeof(struct vector3));
        s->f_bias_eom[i].force   = calloc(1, sizeof(struct vector3));
        s->f_bias_app[i].torque  = calloc(1, sizeof(struct vector3));
        s->f_bias_app[i].force   = calloc(1, sizeof(struct vector3));
        s->f_bias_tf[i].torque   = calloc(1, sizeof(struct vector3));
        s->f_bias_tf[i].force    = calloc(1, sizeof(struct vector3));
        s->f_bias_nact[i].torque = calloc(1, sizeof(struct vector3));
        s->f_bias_nact[i].force  = calloc(1, sizeof(struct vector3));
        // Feed-forward torque
        s->f_ff_app[i].torque = calloc(1, sizeof(struct vector3));
        s->f_ff_app[i].force  = calloc(1, sizeof(struct vector3));
        s->f_ff_jnt[i].torque = calloc(1, sizeof(struct vector3));
        s->f_ff_jnt[i].force  = calloc(1, sizeof(struct vector3));
        // External force
        s->f_ext[i].torque     = calloc(1, sizeof(struct vector3));
        s->f_ext[i].force      = calloc(1, sizeof(struct vector3));
        s->f_ext_app[i].torque = calloc(1, sizeof(struct vector3));
        s->f_ext_app[i].force  = calloc(1, sizeof(struct vector3));
        s->f_ext_tf[i].torque  = calloc(1, sizeof(struct vector3));
        s->f_ext_tf[i].force   = calloc(1, sizeof(struct vector3));
        // Constraint force
        s->f_cstr_app[i].torque = calloc(NR_CSTR, sizeof(struct vector3));
        s->f_cstr_app[i].force  = calloc(NR_CSTR, sizeof(struct vector3));
    }

    for (int i = 0; i < NR_SEGMENTS_WITH_BASE; i++) {
        // FPK
        s->x_tot[i].rotation    = calloc(1, sizeof(struct matrix3x3));
        s->x_tot[i].translation = calloc(1, sizeof(struct vector3));
        // FVK
        s->xd[i].angular_velocity        = calloc(1, sizeof(struct vector3));
        s->xd[i].linear_velocity         = calloc(1, sizeof(struct vector3));
        // FAK
        s->xdd[i].angular_acceleration    = calloc(1, sizeof(struct vector3));
        s->xdd[i].linear_acceleration     = calloc(1, sizeof(struct vector3));
        s->xdd_tf[i].angular_acceleration = calloc(1, sizeof(struct vector3));
        s->xdd_tf[i].linear_acceleration  = calloc(1, sizeof(struct vector3));
        // Inertial force
        s->f_bias_art[i].torque = calloc(1, sizeof(struct vector3));
        s->f_bias_art[i].force  = calloc(1, sizeof(struct vector3));
        // Feed-forward torque
        s->f_ff_art[i].torque = calloc(1, sizeof(struct vector3));
        s->f_ff_art[i].force  = calloc(1, sizeof(struct vector3));
        // External force
        s->f_ext_art[i].torque = calloc(1, sizeof(struct vector3));
        s->f_ext_art[i].force  = calloc(1, sizeof(struct vector3));
        // Constraint force
        s->f_cstr_art[i].torque = calloc(NR_CSTR, sizeof(struct vector3));
        s->f_cstr_art[i].force  = calloc(NR_CSTR, sizeof(struct vector3));
        s->d_cstr_art[i]        = calloc(NR_CSTR * NR_CSTR, sizeof(double));
    }

    for (int i = 0; i < s->nd; i++) {
        // Constraint force
        s->tau_cstr_art[i] = calloc(NR_CSTR, sizeof(joint_torque));
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
    // FVK
    s->xd_jnt = calloc(NR_SEGMENTS, sizeof(struct ga_twist));
    s->xd_tf  = calloc(NR_SEGMENTS, sizeof(struct ga_twist));
    s->xd     = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ga_twist));
    // FAK
    s->xdd_jnt  = calloc(NR_SEGMENTS, sizeof(struct ga_acc_twist));
    s->xdd_net  = calloc(NR_SEGMENTS, sizeof(struct ga_acc_twist));
    s->xdd_bias = calloc(NR_SEGMENTS, sizeof(struct ga_acc_twist));
    s->xdd_tf   = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ga_acc_twist));
    s->xdd_nact = calloc(NR_SEGMENTS, sizeof(struct ga_acc_twist));
    s->xdd      = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ga_acc_twist));
    // Inertia
    s->m_art = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_abi));
    s->m_app = calloc(NR_SEGMENTS, sizeof(struct ma_abi));
    s->m_tf  = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_abi));
    // Inertial force
    s->p           = calloc(NR_SEGMENTS, sizeof(struct ma_momentum));
    s->f_bias_art  = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_wrench));
    s->f_bias_eom  = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_bias_app  = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_bias_tf   = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_bias_nact = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    // Feed-forward torque
    s->f_ff_art = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_wrench));
    s->f_ff_app = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_ff_jnt = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    // External force
    s->f_ext     = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_ext_art = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_wrench));
    s->f_ext_app = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    s->f_ext_tf  = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));
    // Constraint force
    s->f_cstr_art = calloc(NR_SEGMENTS_WITH_BASE, sizeof(struct ma_wrench));
    s->f_cstr_app = calloc(NR_SEGMENTS, sizeof(struct ma_wrench));

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
    // FVK
    s->xd[0].target_body = body_world;
    s->xd[0].reference_body = body_world;
    s->xd[0].point = point_world_origin;
    s->xd[0].frame = frame_world;
    // FAK
    s->xdd[0].target_body = body_world;
    s->xdd[0].reference_body = body_world;
    s->xdd[0].point = point_world_origin;
    s->xdd[0].frame = frame_world;
    // Inertia
    s->m_art[0].body = body_world;
    s->m_art[0].point = point_world_origin;
    s->m_art[0].frame = frame_world;
    // Inertial force
    s->f_bias_art[0].body = body_world;
    s->f_bias_art[0].point = point_world_origin;
    s->f_bias_art[0].frame = frame_world;
    // Feed-forward torque
    s->f_ff_art[EE].body = body_ee;
    s->f_ff_art[EE].point = point_ee_origin;
    s->f_ff_art[EE].frame = frame_ee;
    // External force
    s->f_ext_art[0].body = body_world;
    s->f_ext_art[0].point = point_world_origin;
    s->f_ext_art[0].frame = frame_world;
    // Constraint force
    s->f_cstr_art[EE].body = body_ee;
    s->f_cstr_art[EE].point = point_ee_origin;
    s->f_cstr_art[EE].frame = frame_ee;

    for (int i = 0; i < NR_SEGMENTS; i++) {
        s->f_ext[i].body = kc->segment[i].joint.target_body;
        s->f_ext[i].frame = kc->segment[i].link.root_frame;
        s->f_ext[i].point = kc->segment[i].link.root_frame->origin;
    }
}
