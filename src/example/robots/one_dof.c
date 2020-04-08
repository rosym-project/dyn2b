#include <dyn2b/example/robots.h>


struct one_dof
{
    struct frame joint_0_frame;
    struct point link_0_root_origin;
    struct point link_0_tip_origin;
    struct point link_1_root_origin;
    struct point link_1_tip_origin;
    struct frame link_0_root;
    struct frame link_0_tip;
    struct frame link_1_root;
    struct frame link_1_tip;
    struct body link_0;
    struct body link_1;
};

static struct one_dof one_dof = {
    .joint_0_frame = { .name = "joint_0_frame" },
    .link_0_root_origin = { .name = "link_0_root_origin" },
    .link_0_tip_origin = { .name = "link_0_tip_origin" },
    .link_1_root_origin = { .name = "link_1_root_origin" },
    .link_1_tip_origin = { .name = "link_1_tip_origin" },
    .link_0_root = { .origin = &one_dof.link_0_root_origin, .name = "link_0_root" },
    .link_0_tip = { .origin = &one_dof.link_0_tip_origin, .name = "link_0_tip" },
    .link_1_root = { .origin = &one_dof.link_1_root_origin, .name = "link_1_root" },
    .link_1_tip = { .origin = &one_dof.link_1_tip_origin, .name = "link_1_tip" },
    .link_0 = { .name = "Link0" },
    .link_1 = { .name = "Link1" }
};

struct kca_segment one_dof_robot_segments_a[] = {
    {
        .joint_attachment = {
            .target_frame = &one_dof.link_0_tip,
            .target_body = &one_dof.link_0,
            .reference_frame = &one_dof.link_0_root,
            .reference_body = &one_dof.link_0
        },
        .joint = {
            .target_frame = &one_dof.link_1_root,
            .target_body = &one_dof.link_1,
            .reference_frame = &one_dof.link_0_tip,
            .reference_body = &one_dof.link_0
        },
        .link = {
            .root_frame = &one_dof.link_1_root
        }
    }
};

struct kcc_segment one_dof_robot_segments_c[] = {
    {
        .joint_attachment = {
            .rotation = (struct matrix3x3 [1]) { {
                .row_x = { 1.0, 0.0, 0.0 },
                .row_y = { 0.0, 1.0, 0.0 },
                .row_z = { 0.0, 0.0, 1.0 } }
            },
            .translation = (struct vector3 [1]) { 0.0, 0.0, 0.0 }
        },
        .joint = {
            .type = JOINT_TYPE_REVOLUTE,
            .revolute_joint = {
                .axis = JOINT_AXIS_Z
            }
        },
        .link = {}
    }
};

struct kca_kinematic_chain one_dof_robot_a = {
    .number_of_segments = 1,
    .segment = (struct kca_segment *)one_dof_robot_segments_a
};

struct kcc_kinematic_chain one_dof_robot_c = {
    .number_of_segments = 1,
    .segment = (struct kcc_segment *)one_dof_robot_segments_c
};