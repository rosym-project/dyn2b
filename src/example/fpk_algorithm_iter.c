#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/example/solver_state.h>
#include <dyn2b/example/robots.h>
#include <dyn2b/example/chain_iterator.h>


void fpk_c()
{
    struct kcc_kinematic_chain *kc = &two_dof_robot_c;
    struct solver_state_c s;

    setup_simple_state_c(kc, &s);

    int curr;
    int next;
    struct kcc_joint joint;
    enum joint_type joint_type;
    struct gc_pose x_link;
    bool has_next;
    struct kcc_iterator_nbx iter = {
        .chain = kc,
        .current_index = &curr,
        .next_index = &next,
        .joint = &joint,
        .joint_type = &joint_type,
        .x_link = &x_link,
        .has_next = &has_next
    };
    kcc_iterator_reset(&iter);

    s.q[0] = 1.0;
    s.q[1] = 1.0;


    while (has_next) {
        // Position
        //

        // X_{J,i}
        kcc_joint[joint_type].fpk(&joint, &s.q[curr], &s.x_jnt[curr]);

        // i^X_{i-1} = X_{J,i} X_{T,i}
        gc_pose_compose(&x_link, &s.x_jnt[curr], &s.x_rel[curr]);

        // i^X_0 = i^X_{i-1} {i-1}^X_0
        gc_pose_compose(&s.x_tot[curr], &s.x_rel[curr], &s.x_tot[next]);

        kcc_iterator_next(&iter);
    }

    gc_pose_log(&s.x_tot[s.nbody]);
}


int main(int argc, char **argv)
{
    fpk_c();

    return 0;
}
