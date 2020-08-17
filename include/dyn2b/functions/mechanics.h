#ifndef DYN2B_FUNCTIONS_MECHANICS_H
#define DYN2B_FUNCTIONS_MECHANICS_H

#include <dyn2b/types/mechanics.h>
#include <dyn2b/types/geometry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Add two <m> x <n> arrays of acceleration energy (coordinates).
 *
 * E[i,j = E_1[i,j] + E_2[i,j]
 */
void mc_eacc_add(
        const mc_eacc *e1,
        const mc_eacc *e2,
        mc_eacc *e,
        int m,
        int n);

/**
 * Balance acceleration energy (coordinates).
 *
 * \param decomp The LDL^T decomposition of a <count> x <count> matrix of
 * acceleration energy.
 *
 * \param e <count> x 1 vector of acceleration energy.
 *
 * scale = E_1^{-1} E_2
 */
void mc_eacc_balance(
        const double *decomp,
        const mc_eacc *e,
        double *scale,
        int count);

/**
 * Spatial cross product's dual (coordinates).
 *
 * Xd x* p
 */
void mc_momentum_derive(
        const struct gc_twist *xd,
        const struct mc_momentum *p,
        struct mc_wrench *r);

/**
 * Spatial cross product's dual (ADT).
 *
 * Xd x* p
 */
void ma_momentum_derive(
        const struct ga_twist *xd,
        const struct ma_momentum *p,
        struct ma_wrench *r);


/**
 * Transform an array of <count> wrenches from the pose's target frame to the
 * pose's reference frame (coordinates).
 *
 * X^T F[i]
 */
void mc_wrench_tf_tgt_to_ref(
        const struct gc_pose *x,
        const struct mc_wrench *f,
        struct mc_wrench *r,
        int count);

/**
 * Transform a wrench from the pose's target frame to the pose's reference frame
 * (ADT).
 *
 * X^T F
 */
void ma_wrench_tf_tgt_to_ref(
        const struct ga_pose *x,
        const struct ma_wrench *f,
        struct ma_wrench *r);

/**
 * Invert an array of <count> wrenches (coordinates).
 *
 * -F
 */
void mc_wrench_invert(
        const struct mc_wrench *f,
        struct mc_wrench *r,
        int count);

/**
 * Invert a wrench (ADT).
 *
 * -F
 */
void ma_wrench_invert(
        const struct ma_wrench *f,
        struct ma_wrench *r);

/**
 * Add two arrays of <count> wrenches component-wise (coordinates).
 *
 * F_1[i] + F_2[i]
 */
void mc_wrench_add(
        const struct mc_wrench *f1,
        const struct mc_wrench *f2,
        struct mc_wrench *r,
        int count);

/**
 * Add two wrenches (ADT).
 *
 * F_1 + F_2
 */
void ma_wrench_add(
        const struct ma_wrench *f1,
        const struct ma_wrench *f2,
        struct ma_wrench *r);

/**
 * Substract two arrays of <count> wrenches component-wise (coordinates).
 *
 * F_1[i] - F_2[i]
 */
void mc_wrench_sub(
        const struct mc_wrench *f1,
        const struct mc_wrench *f2,
        struct mc_wrench *r,
        int count);

/**
 * Substract two wrenches (ADT).
 *
 * F_1 - F_2
 */
void ma_wrench_sub(
        const struct ma_wrench *f1,
        const struct ma_wrench *f2,
        struct ma_wrench *r);

/**
 * Print an array of <count> wrenches (coordinates).
 */
void mc_wrench_log(
        const struct mc_wrench *f,
        int count);

/**
 * Print a wrench (ADT).
 */
void ma_wrench_log(
        const struct ma_wrench *f);


/**
 * Map a twist into a momentum (coordinates).
 *
 * M Xd
 */
void mc_rbi_map_twist_to_momentum(
        const struct mc_rbi *m,
        const struct gc_twist *xd,
        struct mc_momentum *r);

/**
 * Map a twist into a momentum (ADT).
 *
 * M Xd
 */
void ma_rbi_map_twist_to_momentum(
        const struct ma_rbi *m,
        const struct ga_twist *xd,
        struct ma_momentum *r);

/**
 * Convert rigid-body inertia to articulated-body inertia (coordinates).
 * 
 * M^A = M
 */
void mc_rbi_to_abi(
        const struct mc_rbi *rbi,
        struct mc_abi *r);

/**
 * Convert rigid-body inertia to articulated-body inertia (ADT).
 * 
 * M^A = M
 */
void ma_rbi_to_abi(
        const struct ma_rbi *rbi,
        struct ma_abi *r);

/**
 * Print rigid-body inertia (coordinates).
 */
void mc_rbi_log(
        const struct mc_rbi *m);

/**
 * Print rigid-body inertia (ADT).
 */
void ma_rbi_log(
        const struct ma_rbi *m);


/**
 * Transform inertia from pose's target frame to pose's reference frame
 * (coordinates).
 * 
 * X^T M^A X
 */
void mc_abi_tf_tgt_to_ref(
        const struct gc_pose *x,
        const struct mc_abi *m,
        struct mc_abi *r);

/**
 * Transform inertia from pose's target frame to pose's reference frame (ADT).
 * 
 * X^T M^A X
 */
void ma_abi_tf_tgt_to_ref(
        const struct ga_pose *x,
        const struct ma_abi *m,
        struct ma_abi *r);

/**
 * Add two articulated-body inertias (coordinates).
 * 
 * M^A_1 + M^A_2
 */
void mc_abi_add(
        const struct mc_abi *m1,
        const struct mc_abi *m2,
        struct mc_abi *r);

/**
 * Add two articulated-body inertias (ADT).
 * 
 * M^A_1 + M^A_2
 */
void ma_abi_add(
        const struct ma_abi *m1,
        const struct ma_abi *m2,
        struct ma_abi *r);

/**
 * Map an acceleration twist into a wrench (coordinates).
 *
 * M^A Xdd
 */
void mc_abi_map_acc_twist_to_wrench(
        const struct mc_abi *m,
        const struct gc_acc_twist *xdd,
        struct mc_wrench *f);

/**
 * Map an acceleration twist into a wrench (ADT).
 *
 * M^A Xdd
 */
void ma_abi_map_acc_twist_to_wrench(
        const struct ma_abi *m,
        const struct ga_acc_twist *xdd,
        struct ma_wrench *f);

/**
 * Print articulated-body inertia (coordinates).
 */
void mc_abi_log(
        const struct mc_abi *m);

/**
 * Print articulated-body inertia (ADT).
 */
void ma_abi_log(
        const struct ma_abi *m);

#ifdef __cplusplus
}
#endif

#endif
