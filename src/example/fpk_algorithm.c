#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/example/solver_state.h>
#include <dyn2b/example/robots.h>


void fpk_a()
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

        // i^X_0 = i^X_{i-1} {i-1}^X_0
        ga_pose_compose(&s.x_rel[i - 1], &s.x_tot[i - 1], &s.x_tot[i]);
    }

    ga_pose_log(&s.x_tot[s.nbody]);
}


void fpk_c()
{
    struct kcc_kinematic_chain *kc = &two_dof_robot_c;
    struct solver_state_c s;

    setup_simple_state_c(kc, &s);

    s.q[0] = 1.0;
    s.q[1] = 1.0;


    for (int i = 1; i < s.nbody + 1; i++) {
        struct kcc_joint *joint = &kc->segment[i - 1].joint;
        int joint_type = joint->type;

        // Position
        //

        // X_{J,i}
        kcc_joint[joint_type].fpk(joint, &s.q[i - 1], &s.x_jnt[i - 1]);

        // i^X_{i-1} = X_{J,i} X_{T,i}
        gc_pose_compose(&s.x_jnt[i - 1], &kc->segment[i - 1].joint_attachment, &s.x_rel[i - 1]);

        // i^X_0 = i^X_{i-1} {i-1}^X_0
        gc_pose_compose(&s.x_rel[i - 1], &s.x_tot[i - 1], &s.x_tot[i]);
    }

    gc_pose_log(&s.x_tot[s.nbody]);
}


int main(int argc, char **argv)
{
    fpk_a();
    fpk_c();

    return 0;
}
