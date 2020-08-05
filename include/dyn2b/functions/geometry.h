#ifndef DYN2B_FUNCTIONS_GEOMETRY_H
#define DYN2B_FUNCTIONS_GEOMETRY_H

#include <dyn2b/types/geometry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Compose two poses (coordinates).
 *
 * X_1 X_2
 */
void gc_pose_compose(
        const struct gc_pose *x1,
        const struct gc_pose *x2,
        struct gc_pose *r);

/**
 * Compose two poses (ADT).
 *
 * X_1 X_2
 */
void ga_pose_compose(
        const struct ga_pose *x1,
        const struct ga_pose *x2,
        struct ga_pose *r);

/**
 * Print a pose (coordinates).
 */
void gc_pose_log(
        const struct gc_pose *x);

/**
 * Print a pose (ADT).
 */
void ga_pose_log(
        const struct ga_pose *x);


/**
 * Transform twist from the pose's reference frame to the pose's target frame
 * (coordinates).
 *
 * X Xd
 */
void gc_twist_tf_ref_to_tgt(
        const struct gc_pose *x,
        const struct gc_twist *xd,
        struct gc_twist *r);

/**
 * Transform twist from the pose's reference frame to the pose's target frame
 * (ADT).
 *
 * X Xd
 */
void ga_twist_tf_ref_to_tgt(
        const struct ga_pose *x,
        const struct ga_twist *xd,
        struct ga_twist *r);

/**
 * Accumulate two twists (coordinates).
 * 
 * Xd_1 + Xd_2
 */
void gc_twist_accumulate(
        const struct gc_twist *xd1,
        const struct gc_twist *xd2,
        struct gc_twist *r);

/**
 * Accumulate two twists (ADT).
 *
 * Xd_1 + Xd_2
 */
void ga_twist_accumulate(
        const struct ga_twist *xd1,
        const struct ga_twist *xd2,
        struct ga_twist *r);

/**
 * Spatial cross product (coordinates).
 *
 * Xd_1 x Xd_2
 */
void gc_twist_derive(
        const struct gc_twist *xd1,
        const struct gc_twist *xd2,
        struct gc_acc_twist *r);

/**
 * Spatial cross product (ADT).
 *
 * Xd_1 x Xd_2
 */
void ga_twist_derive(
        const struct ga_twist *xd1,
        const struct ga_twist *xd2,
        struct ga_acc_twist *r);

/**
 * Print a twist (coordinates).
 */
void gc_twist_log(
        const struct gc_twist *xd);

/**
 * Print a twist (ADT).
 */
void ga_twist_log(
        const struct ga_twist *xd);


/**
 * Transform acceleration twist from the pose's reference frame to the pose's
 * target frame (coordinates).
 *
 * X Xdd
 */
void gc_acc_twist_tf_ref_to_tgt(
        const struct gc_pose *x,
        const struct gc_acc_twist *xdd,
        struct gc_acc_twist *r);

/**
 * Transform acceleration twist from the pose's reference frame to the pose's
 * target frame (ADT).
 *
 * X Xdd
 */
void ga_acc_twist_tf_ref_to_tgt(
        const struct ga_pose *x,
        const struct ga_acc_twist *xdd,
        struct ga_acc_twist *r);

/**
 * Add two acceleration twists (ADT).
 *
 * Xdd_1 + Xdd_2
 */
void ga_acc_twist_add(
        const struct ga_acc_twist *xdd1,
        const struct ga_acc_twist *xdd2,
        struct ga_acc_twist *r);

/**
 * Add two acceleration twists (coordinates).
 *
 * Xdd_1 + Xdd_2
 */
void gc_acc_twist_add(
        const struct gc_acc_twist *xdd1,
        const struct gc_acc_twist *xdd2,
        struct gc_acc_twist *r);

/**
 * Accumulate two acceleration twists (coordinates).
 *
 * Xdd_1 + Xdd_2
 */
void gc_acc_twist_accumulate(
        const struct gc_acc_twist *xdd1,
        const struct gc_acc_twist *xdd2,
        struct gc_acc_twist *r);

/**
 * Accumulate two acceleration twists (ADT).
 *
 * Xdd_1 + Xdd_2
 */
void ga_acc_twist_accumulate(
        const struct ga_acc_twist *xdd1,
        const struct ga_acc_twist *xdd2,
        struct ga_acc_twist *r);

/**
 * Print an acceleration twist (coordinates).
 */
void gc_acc_twist_log(
        const struct gc_acc_twist *xdd);

/**
 * Print an acceleration twist (ADT).
 */
void ga_acc_twist_log(
        const struct ga_acc_twist *xdd);

#ifdef __cplusplus
}
#endif

#endif
