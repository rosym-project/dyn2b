#ifndef DYN2B_TYPES_GEOMETRY_NBX_H
#define DYN2B_TYPES_GEOMETRY_NBX_H

#include <dyn2b/types/geometry.h>

#ifdef __cplusplus
extern "C" {
#endif


struct gc_pose_compose_nbx {
    const struct gc_pose *x1;
    const struct gc_pose *x2;
    struct gc_pose *r;
};

#ifdef __cplusplus
}
#endif

#endif
