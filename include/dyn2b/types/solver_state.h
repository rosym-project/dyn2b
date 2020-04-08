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

    // joint motion state
    joint_position *q;              // joint position                           [nq]
    joint_velocity *qd;             // joint velocity                           [nd]
};

#ifdef __cplusplus
}
#endif

#endif