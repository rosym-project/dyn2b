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

#ifdef __cplusplus
}
#endif

#endif
