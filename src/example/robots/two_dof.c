#include <dyn2b/example/robots.h>


struct two_dof
{
    struct frame joint_0_frame;
    struct frame joint_1_frame;
    struct point link_0_root_origin;
    struct point link_0_tip_origin;
    struct point link_1_root_origin;
    struct point link_1_tip_origin;
    struct point link_2_root_origin;
    struct point link_2_tip_origin;
    struct frame link_0_root;
    struct frame link_0_tip;
    struct frame link_1_root;
    struct frame link_1_tip;
    struct frame link_2_root;
    struct frame link_2_tip;
    struct body link_0;
    struct body link_1;
    struct body link_2;
};

static struct two_dof two_dof = {
    .joint_0_frame = { .name = "joint_0_frame" },
    .joint_1_frame = { .name = "joint_1_frame" },
    .link_0_root_origin = { .name = "link_0_root_origin" },
    .link_0_tip_origin = { .name = "link_0_tip_origin" },
    .link_1_root_origin = { .name = "link_1_root_origin" },
    .link_1_tip_origin = { .name = "link_1_tip_origin" },
    .link_2_root_origin = { .name = "link_2_root_origin" },
    .link_2_tip_origin = { .name = "link_2_tip_origin" },
    .link_0_root = { .origin = &two_dof.link_0_root_origin, .name = "link_0_root" },
    .link_0_tip = { .origin = &two_dof.link_0_tip_origin, .name = "link_0_tip" },
    .link_1_root = { .origin = &two_dof.link_1_root_origin, .name = "link_1_root" },
    .link_1_tip = { .origin = &two_dof.link_1_tip_origin, .name = "link_1_tip" },
    .link_2_root = { .origin = &two_dof.link_2_root_origin, .name = "link_2_root" },
    .link_2_tip = { .origin = &two_dof.link_2_tip_origin, .name = "link_2_tip" },
    .link_0 = { .name = "Link0" },
    .link_1 = { .name = "Link1" },
    .link_2 = { .name = "Link2" }
};

struct kca_segment two_dof_robot_segments_a[] = {
    {
        .joint_attachment = {
            .target_frame = &two_dof.link_0_tip,
            .target_body = &two_dof.link_0,
            .reference_frame = &two_dof.link_0_root,
            .reference_body = &two_dof.link_0
        },
        .joint = {
            .target_frame = &two_dof.link_1_root,
            .target_body = &two_dof.link_1,
            .reference_frame = &two_dof.link_0_tip,
            .reference_body = &two_dof.link_0
        },
        .link = {
            .root_frame = &two_dof.link_1_root,
            .inertia = {
                .frame = &two_dof.link_1_root,
                .point = &two_dof.link_1_root_origin,
                .body = &two_dof.link_1
            }
        }
    },
    {
        .joint_attachment = {
            .target_frame = &two_dof.link_1_tip,
            .target_body = &two_dof.link_1,
            .reference_frame = &two_dof.link_1_root,
            .reference_body = &two_dof.link_1
        },
        .joint = {
            .target_frame = &two_dof.link_2_root,
            .target_body = &two_dof.link_2,
            .reference_frame = &two_dof.link_1_tip,
            .reference_body = &two_dof.link_1
        },
        .link = {
            .root_frame = &two_dof.link_2_root,
            .inertia = {
                .frame = &two_dof.link_2_root,
                .point = &two_dof.link_2_root_origin,
                .body = &two_dof.link_2
            }
        }
    }
};

struct kcc_segment two_dof_robot_segments_c[] = {
    {
        .joint_attachment = {
            .rotation = (struct matrix3x3 [1]) { {
                .row_x = { 1.0, 0.0, 0.0 },
                .row_y = { 0.0, 1.0, 0.0 },
                .row_z = { 0.0, 0.0, 1.0 }
            } },
            .translation = (struct vector3 [1]) { 0.0, 0.0, 0.0 }
        },
        .joint = {
            .type = JOINT_TYPE_REVOLUTE,
            .revolute_joint = {
                .axis = JOINT_AXIS_Z,
                .inertia = (double[]) { 1.0 }
            }
        },
        .link = {
            .inertia = {
                .zeroth_moment_of_mass = 2.0,
                .first_moment_of_mass = { 2.0, 0.0, 0.0 },      // [1.0, 0.0, 0.0] * 2.0
                .second_moment_of_mass = {
                    .row_x = { 0.0, 0.0, 0.0 },
                    .row_y = { 0.0, 2.0, 0.0 },
                    .row_z = { 0.0, 0.0, 2.0 }
                }
            }
        }
    },
    {
        .joint_attachment = {
            .rotation = (struct matrix3x3 [1]) { {
                .row_x = { 1.0, 0.0, 0.0 },
                .row_y = { 0.0, 1.0, 0.0 },
                .row_z = { 0.0, 0.0, 1.0 }
            } },
            .translation = (struct vector3 [1]) { 2.0, 0.0, 0.0 }
        },
        .joint = {
            .type = JOINT_TYPE_REVOLUTE,
            .revolute_joint = {
                .axis = JOINT_AXIS_Z,
                .inertia = (double[]) { 1.0 }
            }
        },
        .link = {
            .inertia = {
                .zeroth_moment_of_mass = 2.0,
                .first_moment_of_mass = { 2.0, 0.0, 0.0 },      // [1.0, 0.0, 0.0] * 2.0
                .second_moment_of_mass = {
                    .row_x = { 0.0, 0.0, 0.0 },
                    .row_y = { 0.0, 2.0, 0.0 },
                    .row_z = { 0.0, 0.0, 2.0 }
                }
            }
        }
    }
};

struct kca_kinematic_chain two_dof_robot_a = {
    .number_of_segments = 2,
    .segment = (struct kca_segment *)two_dof_robot_segments_a
};

struct kcc_kinematic_chain two_dof_robot_c = {
    .number_of_segments = 2,
    .segment = (struct kcc_segment *)two_dof_robot_segments_c
};