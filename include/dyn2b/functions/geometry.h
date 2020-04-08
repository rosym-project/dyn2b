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
 * Print a twist (coordinates).
 */
void gc_twist_log(
        const struct gc_twist *xd);

/**
 * Print a twist (ADT).
 */
void ga_twist_log(
        const struct ga_twist *xd);

#ifdef __cplusplus
}
#endif

#endif
