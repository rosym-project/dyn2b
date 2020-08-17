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
 * kca_*: _a_bstract data type i.e. the coordinate-independent representation
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

/**
 * Inverse force kinematics (ADT).
 *
 * tau = S^T F
 */
void kca_ifk(
        const struct kca_joint *joint,
        const struct ma_wrench *f);

/**
 * Forward force dynamics (ADT).
 *
 * F = M S D^{-1} tau
 */
void kca_ffd(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        struct ma_wrench *f);

/**
 * Project inertia over a joint (ADT).
 *
 * M^a = P^T M^A
 */
void kca_project_inertia(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        struct ma_abi *r);

/**
 * Project wrench over a joint (ADT).
 *
 * F^a = P^T F^A
 */
void kca_project_wrench(
        const struct kca_joint *joint,
        const struct ma_abi *m,
        const struct ma_wrench *f,
        struct ma_wrench *r);


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

    /**
     * Inverse force kinematics on <count> wrenches (coordinates).
     *
     * tau[i] = S^T F[i]
     */
    void (*ifk)(
            const struct kcc_joint *joint,
            const struct mc_wrench *f,
            joint_torque *tau,
            int count);

    /**
     * Forward force dynamics on <count> torques (coordinates).
     *
     * F[i] = M S D^{-1} tau[i]
     */
    void (*ffd)(
            const struct kcc_joint *joint,
            const struct mc_abi *m,
            const joint_torque *tau,
            struct mc_wrench *f,
            int count);

    /**
     * Project inertia over a joint (coordinates).
     *
     * M^a = P^T M^A = (1 - M^A S D^{-1} S^T) M^A
     */
    void (*project_inertia)(
            const struct kcc_joint *joint,
            const struct mc_abi *m,
            struct mc_abi *r);

    /**
     * Project an array of <count> wrenches over a joint (coordinates).
     *
     * F^a[i] = P^T F^A[i] = (1 - M^A S D^{-1} S^T) F^A[i]
     */
    void (*project_wrench)(
            const struct kcc_joint *joint,
            const struct mc_abi *m,
            const struct mc_wrench *f,
            struct mc_wrench *r,
            int count);

    /**
     * Decompose auto constraint energy.
     *
     * E_cstr^A + (F_cstr^A)^T S D^{-1} S^T F_cstr^A
     */
    void (*decomp_e_cstr)(
            const struct kcc_joint *joint,
            const struct mc_abi *m,
            const struct mc_wrench *f,
            double *d_cstr,
            double *r,
            int count);
};


extern const struct kcc_joint_operators kcc_joint[];

#ifdef __cplusplus
}
#endif

#endif
