#ifndef DYN2B_EXAMPLE_ROBOTS_H
#define DYN2B_EXAMPLE_ROBOTS_H

#include <dyn2b/types/kinematic_chain.h>

#ifdef __cplusplus
extern "C" {
#endif


struct kcc_kinematic_chain one_dof_robot_c;
struct kca_kinematic_chain one_dof_robot_a;

struct kcc_kinematic_chain two_dof_robot_c;
struct kca_kinematic_chain two_dof_robot_a;

#ifdef __cplusplus
}
#endif

#endif
