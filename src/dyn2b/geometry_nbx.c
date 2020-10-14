#include <dyn2b/functions/geometry_nbx.h>
#include <dyn2b/functions/geometry.h>

#include <assert.h>


void gc_pose_compose_nbx(
        struct gc_pose_compose_nbx *nbx)
{
    assert(nbx);

    gc_pose_compose(nbx->x1, nbx->x2, nbx->r);
}

