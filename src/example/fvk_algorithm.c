#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/example/solver_state.h>
#include <dyn2b/example/robots.h>


void fvk_a()
{
    struct kca_kinematic_chain *kc = &two_dof_robot_a;
    struct solver_state_a s;

    setup_simple_state_a(kc, &s);


    for (int i = 1; i < s.nbody + 1; i++) {
        struct kca_joint *joint = &kc->segment[i - 1].joint;

        // Position
        //

        // X_{J,i}
        kca_fpk(joint, &s.x_jnt[i - 1]);

        // i^X_{i-1} = X_{J,i} X_{T,i}
        ga_pose_compose(&s.x_jnt[i - 1], &kc->segment[i - 1].joint_attachment, &s.x_rel[i - 1]);


        // Velocity
        //

        // Xd_{J,i} = S qd
        kca_fvk(joint, &s.xd_jnt[i - 1]);

        // Xd_{i-1}' = i^X_{i-1} Xd_{i-1}
        ga_twist_tf_ref_to_tgt(&s.x_rel[i - 1], &s.xd[i - 1], &s.xd_tf[i - 1]);

        // Xd_i = Xd_{i-1}' + Xd_{J,i}
        ga_twist_accumulate(&s.xd_tf[i - 1], &s.xd_jnt[i - 1], &s.xd[i]);
    }

    ga_twist_log(&s.xd[s.nbody]);
}


void fvk_c()
{
    struct kcc_kinematic_chain *kc = &two_dof_robot_c;
    struct solver_state_c s;

    setup_simple_state_c(kc, &s);

    s.q[0] = 1.0;
    s.q[1] = 1.0;
    s.qd[0] = 1.0;
    s.qd[1] = 1.0;


    for (int i = 1; i < s.nbody + 1; i++) {
        struct kcc_joint *joint = &kc->segment[i - 1].joint;
        int joint_type = joint->type;

        // Position
        //

        // X_{J,i}
        kcc_joint[joint_type].fpk(joint, &s.q[i - 1], &s.x_jnt[i - 1]);

        // i^X_{i-1} = X_{J,i} X_{T,i}
        gc_pose_compose(&s.x_jnt[i - 1], &kc->segment[i - 1].joint_attachment, &s.x_rel[i - 1]);


        // Velocity
        //

        // Xd_{J,i} = S qd
        kcc_joint[joint_type].fvk(joint, &s.qd[i - 1], &s.xd_jnt[i - 1]);

        // Xd_{i-1}' = i^X_{i-1} Xd_{i-1}
        gc_twist_tf_ref_to_tgt(&s.x_rel[i - 1], &s.xd[i - 1], &s.xd_tf[i - 1]);

        // Xd_i = Xd_{i-1}' + Xd_{J,i}
        gc_twist_accumulate(&s.xd_tf[i - 1], &s.xd_jnt[i - 1], &s.xd[i]);
    }

    gc_twist_log(&s.xd[s.nbody]);
}


int main(int argc, char **argv)
{
    fvk_a();
    fvk_c();

    return 0;
}
