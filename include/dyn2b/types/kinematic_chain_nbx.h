#ifndef DYN2B_TYPES_KINEMATIC_CHAIN_NBX_H
#define DYN2B_TYPES_KINEMATIC_CHAIN_NBX_H

#include <dyn2b/types/kinematic_chain.h>

#ifdef __cplusplus
extern "C" {
#endif

struct kcc_fpk_nbx {
    const struct kcc_joint *joint;
    const joint_position *q;
    struct gc_pose *x;
};

#ifdef __cplusplus
}
#endif

#endif
