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
};


extern const struct kcc_joint_operators kcc_joint[];

#ifdef __cplusplus
}
#endif

#endif
