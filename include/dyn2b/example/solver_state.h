#ifndef DYN2B_EXAMPLE_SOLVER_STATE_H
#define DYN2B_EXAMPLE_SOLVER_STATE_H

#include <dyn2b/types/solver_state.h>
#include <dyn2b/types/kinematic_chain.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Setup the solver state for a simple task and a serial kinematic chain.
 */
void setup_simple_state_c(
        const struct kcc_kinematic_chain *kc,
        struct solver_state_c *s);

/**
 * Setup the solver state for a simple task and a serial kinematic chain.
 */
void setup_simple_state_a(
        const struct kca_kinematic_chain *kc,
        struct solver_state_a *s);

#ifdef __cplusplus
}
#endif

#endif
