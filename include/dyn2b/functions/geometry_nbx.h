#ifndef DYN2B_FUNCTIONS_GEOMETRY_NBX_H
#define DYN2B_FUNCTIONS_GEOMETRY_NBX_H

#include <dyn2b/types/geometry_nbx.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Compose two poses (coordinates).
 *
 * X_1 X_2
 */
void gc_pose_compose_nbx(
        struct gc_pose_compose_nbx *nbx);


#ifdef __cplusplus
}
#endif

#endif
