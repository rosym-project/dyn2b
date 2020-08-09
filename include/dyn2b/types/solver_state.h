#ifndef DYN2B_TYPES_SOLVER_STATE_H
#define DYN2B_TYPES_SOLVER_STATE_H

#include <dyn2b/types/geometry.h>
#include <dyn2b/types/kinematic_chain.h>

#ifdef __cplusplus
extern "C" {
#endif


struct solver_state_a
{
    int nbody;                      // number of bodies
    int nq;                         // number of joint positions
    int nd;                         // number of motion DoFs

    // spatial motion state
    struct ga_pose      *x_jnt;     // pose over the joint                      [nbody]
    struct ga_pose      *x_rel;     // pose w.r.t. to predecessor body          [nbody]
    struct ga_pose      *x_tot;     // pose relative to base                    [nbody]
    struct ga_twist     *xd_jnt;    // velocity over the joint                  [nbody]
    struct ga_twist     *xd_tf;     // transformed predecessor velocity         [nbody]
    struct ga_twist     *xd;        // velocity                                 [nbody]
    struct ga_acc_twist *xdd_jnt;   // acceleration over the joint              [nbody]
    struct ga_acc_twist *xdd_bias;  // bias acceleration                        [nbody]
    struct ga_acc_twist *xdd_net;   // net acceleration over the joint          [nbody]
    struct ga_acc_twist *xdd_tf;    // tf'ed acceleration                       [nbody]
    struct ga_acc_twist *xdd_nact;  // acceleration w/o active joint contrib.   [nbody]
    struct ga_acc_twist *xdd;       // acceleration                             [nbody]

    // inertia
    struct ma_abi *m_art;           // articulated-body inertia                 [nbody]
    struct ma_abi *m_app;           // apparent inertia                         [nbody]
    struct ma_abi *m_tf;            // tf'ed apparent inertia                   [nbody]

    // inertial force
    struct ma_momentum *p;          // momentum                                 [nbody]
    struct ma_wrench *f_bias_art;   // articulated bias force                   [nbody]
    struct ma_wrench *f_bias_eom;   // articulated equation of motion           [nbody]
    struct ma_wrench *f_bias_app;   // apparent bias force                      [nbody]
    struct ma_wrench *f_bias_tf;    // tf'ed apparent bias force                [nbody]
    struct ma_wrench *f_bias_nact;  // inertial force w/o active joint contrib. [nbody]

    // feed-forward joint torque motion driver
    struct ma_wrench *f_ff_art;     // articulated feed-forward force           [nbody]
    struct ma_wrench *f_ff_app;     // apparent feed-forward force              [nbody]
    struct ma_wrench *f_ff_jnt;     // joint contribution to feed-forward force [nbody]

    // external force motion driver
    struct ma_wrench *f_ext;        // external force                           [nbody]
    struct ma_wrench *f_ext_art;    // articulated external force               [nbody]
    struct ma_wrench *f_ext_app;    // apparent external force                  [nbody]
    struct ma_wrench *f_ext_tf;     // tf'ed apparent external force            [nbody]

    // acceleration constraint motion driver
    struct ma_wrench *f_cstr;       // constraint force                         [nc]
    struct ma_wrench *f_cstr_art;   // articulated constraint force             [nc * nbody]
    struct ma_wrench *f_cstr_app;   // apparent constraint force                [nc * nbody]
};


struct solver_state_c
{
    int nbody;                      // number of bodies
    int nq;                         // number of joint positions
    int nd;                         // number of motion DoFs

    // spatial motion state
    struct gc_pose      *x_jnt;     // pose over the joint                      [nbody]
    struct gc_pose      *x_rel;     // pose w.r.t. to predecessor body          [nbody]
    struct gc_pose      *x_tot;     // pose relative to base                    [nbody]
    struct gc_twist     *xd_jnt;    // velocity over the joint                  [nbody]
    struct gc_twist     *xd_tf;     // transformed predecessor velocity         [nbody]
    struct gc_twist     *xd;        // velocity                                 [nbody]
    struct gc_acc_twist *xdd_jnt;   // acceleration over the joint              [nbody]
    struct gc_acc_twist *xdd_bias;  // bias acceleration                        [nbody]
    struct gc_acc_twist *xdd_net;   // net acceleration over the joint          [nbody]
    struct gc_acc_twist *xdd_tf;    // tf'ed acceleration                       [nbody]
    struct gc_acc_twist *xdd_nact;  // acceleration w/o active joint contrib.   [nbody]
    struct gc_acc_twist *xdd;       // acceleration                             [nbody]

    // joint motion state
    joint_position *q;              // joint position                           [nq]
    joint_velocity *qd;             // joint velocity                           [nd]
    joint_acceleration *qdd;        // joint acceleration                       [nd]

    // inertia
    struct mc_abi *m_art;           // articulated-body inertia                 [nbody]
    struct mc_abi *m_app;           // apparent inertia                         [nbody]
    struct mc_abi *m_tf;            // tf'ed apparent inertia                   [nbody]
    joint_inertia *d;               // constrained inertia                      [nd]

    // inertial force
    struct mc_momentum *p;          // momentum                                 [nbody]
    struct mc_wrench *f_bias_art;   // articulated bias force                   [nbody]
    struct mc_wrench *f_bias_eom;   // articulated equation of motion           [nbody]
    struct mc_wrench *f_bias_app;   // apparent bias force                      [nbody]
    struct mc_wrench *f_bias_tf;    // tf'ed apparent bias force                [nbody]
    struct mc_wrench *f_bias_nact;  // inertial force w/o active joint contrib. [nbody]
    joint_torque *tau_bias_art;     // torque due to art. bias force            [nd]

    // feed-forward joint torque motion driver
    joint_torque *tau_ff;           // feed-forward torque                      [nd]
    struct mc_wrench *f_ff_art;     // articulated feed-forward force           [nbody]
    struct mc_wrench *f_ff_app;     // apparent feed-forward force              [nbody]
    struct mc_wrench *f_ff_jnt;     // joint contribution to feed-forward force [nbody]
    joint_torque *tau_ff_art;       // torque due to art. feed-forward force    [nd]

    // external force motion driver
    struct mc_wrench *f_ext;        // external force                           [nbody]
    struct mc_wrench *f_ext_art;    // articulated external force               [nbody]
    struct mc_wrench *f_ext_app;    // apparent external force                  [nbody]
    struct mc_wrench *f_ext_tf;     // tf'ed apparent external force            [nbody]
    joint_torque *tau_ext_art;      // torque due to art. external force        [nd]

    // acceleration constraint motion driver
    struct mc_wrench *f_cstr;       // constraint force                         [nc]
    struct mc_wrench *f_cstr_art;   // articulated constraint force             [nbody * nc]
    struct mc_wrench *f_cstr_app;   // apparent constraint force                [nbody * nc]
    joint_torque **tau_cstr_art;    // torque due to art. cstr. force           [nd * nc]

    joint_torque *tau_ctrl;         // joint control torque                     [nd]
};

#ifdef __cplusplus
}
#endif

#endif
