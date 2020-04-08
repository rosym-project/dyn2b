#ifndef DYN2B_TYPES_KINEMATIC_CHAIN_H
#define DYN2B_TYPES_KINEMATIC_CHAIN_H

#include <dyn2b/types/geometry.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef double joint_position;


enum joint_type
{
    JOINT_TYPE_REVOLUTE
};

enum joint_axis
{
    JOINT_AXIS_X = 0,
    JOINT_AXIS_Y = 1,
    JOINT_AXIS_Z = 2
};

struct kcc_revolute_joint
{
    enum joint_axis axis;
};

struct kcc_joint
{
    enum joint_type type;

    union
    {
        struct kcc_revolute_joint revolute_joint;
    };
};

struct kca_joint
{
    struct frame *target_frame;
    struct body *target_body;
    struct frame *reference_frame;
    struct body *reference_body;
};

struct kcc_link
{
};

struct kca_link
{
    struct frame *root_frame;
};

struct kcc_segment
{
    struct gc_pose joint_attachment;
    struct kcc_joint joint;
    struct kcc_link link;
};

struct kca_segment
{
    struct ga_pose joint_attachment;
    struct kca_joint joint;
    struct kca_link link;
};

struct kca_kinematic_chain
{
    int number_of_segments;
    struct kca_segment *segment;
};

struct kcc_kinematic_chain
{
    int number_of_segments;
    struct kcc_segment *segment;
};

#ifdef __cplusplus
}
#endif

#endif
