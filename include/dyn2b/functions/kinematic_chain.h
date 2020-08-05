#ifndef DYN2B_FUNCTIONS_KINEMATIC_CHAIN_H
#define DYN2B_FUNCTIONS_KINEMATIC_CHAIN_H

#include <dyn2b/types/kinematic_chain.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Data types to represent rigid body systems (kinematic chains).
 * 
 * The following naming conventions apply:
 * kcc_*: _c_oordinate representation
 * kcaa_*: _a_bstract data type i.e. the coordinate-independent representation
 */


/**
 * Forward position kinematics (ADT).
 * 
 * X_J
 */
void kca_fpk(
        const struct kca_joint *joint,
        struct ga_pose *x);

/**
 * Forward velocity kinematics (ADT).
 * 
 * Xd = S qd
 */
void kca_fvk(
        const struct kca_joint *joint,
        struct ga_twist *xd);

/**
 * Foward acceleration kinematics (ADT).
 *
 * Xdd = S qdd
 */
void kca_fak(
        const struct kca_joint *joint,
        struct ga_acc_twist *xdd);

/**
 * Velocity-dependent inertial acceleration of the joint's target body (ADT).
 *
 * Xdd_{bias} = Sd qd + (Xd x S) qd
 */
void kca_inertial_acceleration(
        const struct kca_joint *joint,
        const struct ga_twist *xd,
        struct ga_acc_twist *xdd);


struct kcc_joint_operators
{
    /**
     * Forward position kinematics (coordinates).
     * 
     * X_J
     */
    void (*fpk)(
            const struct kcc_joint *joint,
            const joint_position *q,
            struct gc_pose *x);

    /**
     * Forward velocity kinematics (coordinates).
     * 
     * Xd = S qd
     */
    void (*fvk)(
            const struct kcc_joint *joint,
            const joint_velocity *qd,
            struct gc_twist *xd);

    /**
     * Forward acceleration kinematics (coordinates).
     *
     * Xdd = S qdd
     */
    void (*fak)(
            const struct kcc_joint *joint,
            const joint_acceleration *qdd,
            struct gc_acc_twist *xdd);

    /**
     * Velocity-dependent inertial acceleration of the joint's target body
     * (coordinates).
     *
     * Xdd_{bias} = Sd qd + (Xd x S) qd
     */
    void (*inertial_acceleration)(
            const struct kcc_joint *joint,
            const struct gc_twist *xd,
            const joint_velocity *qd,
            struct gc_acc_twist *xdd);
};


extern const struct kcc_joint_operators kcc_joint[];

#ifdef __cplusplus
}
#endif

#endif
