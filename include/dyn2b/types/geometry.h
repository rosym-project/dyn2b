#ifndef DYN2B_TYPES_GEOMETRY_H
#define DYN2B_TYPES_GEOMETRY_H

#include <dyn2b/types/linear_algebra.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * Data types to represent geometric aspects of rigid bodies.
 * 
 * The following naming conventions apply:
 * gc_*: _c_oordinate representation
 * ga_*: _a_bstract data type i.e. the coordinate-independent representation
 */


struct point
{
    const char *name;
};

struct frame
{
    struct point *origin;
    const char *name;
};

struct body
{
    const char *name;
};


struct gc_pose
{
    struct matrix3x3 *rotation;
    struct vector3 *translation;
};

struct ga_pose
{
    struct body *target_body;
    struct body *reference_body;
    struct frame *target_frame;
    struct frame *reference_frame;
};


struct gc_twist
{
    struct vector3 *angular_velocity;
    struct vector3 *linear_velocity;
};

struct ga_twist
{
    struct body *target_body;
    struct body *reference_body;
    struct point *point;
    struct frame *frame;
};


struct gc_acc_twist
{
    struct vector3 *angular_acceleration;
    struct vector3 *linear_acceleration;
};

struct ga_acc_twist
{
    struct body *target_body;
    struct body *reference_body;
    struct point *point;
    struct frame *frame;
};

#ifdef __cplusplus
}
#endif

#endif
