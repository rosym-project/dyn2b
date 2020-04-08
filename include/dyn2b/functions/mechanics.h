#ifndef DYN2B_FUNCTIONS_MECHANICS_H
#define DYN2B_FUNCTIONS_MECHANICS_H

#include <dyn2b/types/mechanics.h>
#include <dyn2b/types/geometry.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Print rigid-body inertia (coordinates).
 */
void mc_rbi_log(
        const struct mc_rbi *m);

/**
 * Print rigid-body inertia (ADT).
 */
void ma_rbi_log(
        const struct ma_rbi *m);

#ifdef __cplusplus
}
#endif

#endif
