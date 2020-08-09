#ifndef DYN2B_TYPES_MECHANICS_H
#define DYN2B_TYPES_MECHANICS_H

#include <dyn2b/types/linear_algebra.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Data types to represent the dynamics of rigid bodies.
 * 
 * The following naming conventions apply:
 * mc_*: _c_oordinate representation
 * ma_*: _a_bstract data type i.e. the coordinate-independent representation
 */


/**
 * Acceleration energy (coordinates)
 */
typedef double mc_eacc;


/**
 * Momentum (coordinates)
 */
struct mc_momentum
{
    struct vector3 *angular_momentum;
    struct vector3 *linear_momentum;
};

/**
 * Momentum (ADT)
 */
struct ma_momentum
{
    struct body *body;
    struct point *point;
    struct frame *frame;
};


/**
 * Wrench (coordinates)
 */
struct mc_wrench
{
    struct vector3 *torque;
    struct vector3 *force;
};

/**
 * Wrench (ADT)
 */
struct ma_wrench
{
    struct body *body;
    struct point *point;
    struct frame *frame;
};


/**
 * Rigid-body inertia (coordinates)
 */
struct mc_rbi
{
    double zeroth_moment_of_mass;           // mass
    struct vector3 first_moment_of_mass;    // mass * centre of mass
    struct matrix3x3 second_moment_of_mass; // rotational inertia
};

/**
 * Rigid-body inertia (ADT)
 */
struct ma_rbi
{
    struct body *body;
    struct point *point;
    struct frame *frame;
};


/**
 * Articulated-body inertia (coordinates)
 */
struct mc_abi
{
    struct matrix3x3 zeroth_moment_of_mass; // "mass"
    struct matrix3x3 first_moment_of_mass;  // "mass * centre of mass"
    struct matrix3x3 second_moment_of_mass; // "rotational inertia"
};

/**
 * Articulated-body inertia (ADT)
 */
struct ma_abi
{
    struct body *body;
    struct point *point;
    struct frame *frame;
};

#ifdef __cplusplus
}
#endif

#endif
