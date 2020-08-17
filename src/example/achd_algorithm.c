#include <dyn2b/functions/geometry.h>
#include <dyn2b/functions/mechanics.h>
#include <dyn2b/functions/kinematic_chain.h>
#include <dyn2b/example/solver_state.h>
#include <dyn2b/example/robots.h>
#include <stdio.h>


void achd_a()
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


        // Acceleration
        //

        // Xdd_{bias,i} = S_i qdd_i + Sd_i qd_i + Xd_i x S_i qd_i
        kca_inertial_acceleration(joint, &s.xd[i], &s.xdd_bias[i - 1]);


        // Inertia
        //

        // M_i^A = M_i
        ma_rbi_to_abi(&kc->segment[i - 1].link.inertia, &s.m_art[i]);


        // Force
        //

        // p_i = M_i Xd_i
        ma_rbi_map_twist_to_momentum(&kc->segment[i - 1].link.inertia, &s.xd[i], &s.p[i - 1]);

        // F_{bias,i}^A = Xd_i x* P_i
        ma_momentum_derive(&s.xd[i], &s.p[i - 1], &s.f_bias_art[i]);


        // F_{ext,i}^A = -F_{ext,i}
        ma_wrench_invert(&s.f_ext[i - 1], &s.f_ext_art[i]);
    }


    for (int i = s.nbody; i > 0; i--) {
        struct kca_joint *joint = &kc->segment[i - 1].joint;

        // Inertia
        //

        // M_i^a = P_i^T M_i^A
        kca_project_inertia(joint, &s.m_art[i], &s.m_app[i - 1]);

        // M_{i-1}^a' = {i-1}^X_i* M_i^a i^X_{i-1}
        ma_abi_tf_tgt_to_ref(&s.x_rel[i - 1], &s.m_app[i - 1], &s.m_tf[i]);

        // M_{i-1}^A += M_{i-1}^a'
        ma_abi_add(&s.m_art[i - 1], &s.m_tf[i], &s.m_art[i - 1]);


        // Force
        //

        // F_{bias,i}^A' = M_i^A Xdd_{bias,i}
        ma_abi_map_acc_twist_to_wrench(&s.m_art[i], &s.xdd_bias[i - 1], &s.f_bias_eom[i - 1]);

        // F_{bias,i}^A' += F_{bias,i}^A
        ma_wrench_add(&s.f_bias_eom[i - 1], &s.f_bias_art[i], &s.f_bias_eom[i - 1]);

        // F_{bias,i}^a = P_i^T F_{bias,i}^A'
        kca_project_wrench(joint, &s.m_art[i], &s.f_bias_eom[i - 1], &s.f_bias_app[i - 1]);

        // F_{bias,i}^a' = {i-1}^X_i* F_{bias,i}^a
        ma_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_bias_app[i - 1], &s.f_bias_tf[i - 1]);

        // F_{bias,i-1}^A += F_{bias,i}^a'
        ma_wrench_add(&s.f_bias_art[i - 1], &s.f_bias_tf[i - 1], &s.f_bias_art[i - 1]);


        // F_{tau,i} = M_i^A S_i D^{-1} tau_{ff,i}
        kca_ffd(joint, &s.m_art[i], &s.f_ff_jnt[i - 1]);

        // F_{tau,i}^a = P_i^T F_{tau,i}^A
        kca_project_wrench(joint, &s.m_art[i], &s.f_ff_art[i], &s.f_ff_app[i - 1]);

        // F_{tau,i}^a += F_{tau,i}'
        ma_wrench_add(&s.f_ff_app[i - 1], &s.f_ff_jnt[i - 1], &s.f_ff_app[i - 1]);

        // F_{tau,i-1}^A = {i-1}^X_i* F_{tau,i}^a
        ma_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_ff_app[i - 1], &s.f_ff_art[i - 1]);


        // F_{ext,i}^a = P_i^T F_{ext,i}^A
        kca_project_wrench(joint, &s.m_art[i], &s.f_ext_art[i], &s.f_ext_app[i - 1]);

        // F_{ext,i}^a' = {i-1}^X_i* F_{ext,i}^a
        ma_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_ext_app[i - 1], &s.f_ext_tf[i - 1]);

        // F_{ext,i-1}^A += F_{ext,i}^a'
        ma_wrench_add(&s.f_ext_art[i - 1], &s.f_ext_tf[i - 1], &s.f_ext_art[i - 1]);


        // F_{cstr,i}^a = P_i^T
        kca_project_wrench(joint, &s.m_art[i], &s.f_cstr_art[i], &s.f_cstr_app[i - 1]);

        // F_{cstr,i-1}^A = {i-1}^X_i* F_{cstr,i}^a
        ma_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_cstr_app[i - 1], &s.f_cstr_art[i - 1]);
    }

    ma_abi_log(&s.m_art[0]);
    ma_wrench_log(&s.f_bias_art[0]);
    ma_wrench_log(&s.f_ff_art[0]);
    ma_wrench_log(&s.f_ext_art[0]);
    ma_wrench_log(&s.f_cstr_art[0]);


    for (int i = 1; i < s.nbody + 1; i++) {
        struct kca_joint *joint = &kc->segment[i - 1].joint;

        // Acceleration
        //

        // Xdd_{i-1}' = i^X_{i-1} Xdd_{i-1}
        ga_acc_twist_tf_ref_to_tgt(&s.x_rel[i - 1], &s.xdd[i - 1], &s.xdd_tf[i - 1]);

        // Xdd_{nact,i-1} = Xdd_{i-1}' + Xdd_{bias,i-1}
        ga_acc_twist_accumulate(&s.xdd_tf[i - 1], &s.xdd_bias[i - 1], &s.xdd_nact[i - 1]);


        // Force
        //

        // F_{nact,i} = M_i^A Xdd_{nact,i-1}'
        ma_abi_map_acc_twist_to_wrench(&s.m_art[i], &s.xdd_nact[i - 1], &s.f_bias_nact[i - 1]);

        // tau_{bias,i} = S^T F_{bias,i}
        kca_ifk(joint, &s.f_bias_nact[i - 1]);


        // tau_{ff,i}^A = S^T F_{ff,i}
        kca_ifk(joint, &s.f_ff_art[i]);


        // tau_{ext,i}^A = S^T F_{ext,i}
        kca_ifk(joint, &s.f_ext_art[i]);


        // tau_{cstr,i}^A' = S^T F_{cstr,i}
        kca_ifk(joint, &s.f_cstr_art[i]);


        // Resultant acceleration
        //

        // Xdd_{J,i} = S_i qdd_i
        kca_fak(joint, &s.xdd_jnt[i - 1]);

        // Xdd_{J,i}' = Xdd_{J,i} + Xdd_{bias,i}
        ga_acc_twist_add(&s.xdd_jnt[i - 1], &s.xdd_bias[i - 1], &s.xdd_net[i - 1]);

        // Xdd_i = Xdd_{i-1}' + Xdd_{J,i}'
        ga_acc_twist_accumulate(&s.xdd_tf[i - 1], &s.xdd_net[i - 1], &s.xdd[i]);
    }

    ga_acc_twist_log(&s.xdd[s.nbody]);
}


void achd_c()
{
    struct kcc_kinematic_chain *kc = &two_dof_robot_c;
    struct solver_state_c s;

    const int NR_CSTR = 6;

    setup_simple_state_c(kc, &s);

    s.q[0] = 1.0;
    s.q[1] = 1.0;
    s.qd[0] = 1.0;
    s.qd[1] = 1.0;
    s.f_ext[s.nbody - 1].force[0].x = 1.0;
    s.f_ext[s.nbody - 1].force[0].y = 1.0;
    s.f_ext[s.nbody - 1].force[0].z = 1.0;
    s.tau_ff[0] = 1.0;
    s.tau_ff[1] = 1.0;
    s.xdd->linear_acceleration[0].z = 9.81;
    s.f_cstr_art[s.nbody].torque[0].x = 1.0;
    s.f_cstr_art[s.nbody].torque[1].y = 1.0;
    s.f_cstr_art[s.nbody].torque[2].z = 1.0;
    s.f_cstr_art[s.nbody].force[3].x = 1.0;
    s.f_cstr_art[s.nbody].force[4].y = 1.0;
    s.f_cstr_art[s.nbody].force[5].z = 1.0;
    s.e_cstr[0] = 1.0;
    s.e_cstr[1] = 1.0;
    s.e_cstr[2] = 1.0;
    s.e_cstr[3] = 1.0;
    s.e_cstr[4] = 1.0;
    s.e_cstr[5] = 1.0;


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


        // Acceleration
        //

        // Xdd_{bias,i} = S_i qdd_i + Sd_i qd_i + Xd_i x S_i qd_i
        kcc_joint[joint_type].inertial_acceleration(joint, &s.xd[i], &s.qd[i - 1], &s.xdd_bias[i - 1]);


        // Inertia
        //

        // M_i^A = M_i
        mc_rbi_to_abi(&kc->segment[i - 1].link.inertia, &s.m_art[i]);


        // Force
        //

        // P_i = M_i Xd_i
        mc_rbi_map_twist_to_momentum(&kc->segment[i - 1].link.inertia, &s.xd[i], &s.p[i - 1]);

        // F_{bias,i}^A = Xd_i x* P_i
        mc_momentum_derive(&s.xd[i], &s.p[i - 1], &s.f_bias_art[i]);


        // F_{ext,i}^A = -F_{ext,i}
        mc_wrench_invert(&s.f_ext[i - 1], &s.f_ext_art[i], 1);
    }


    for (int i = s.nbody; i > 0; i--) {
        struct kcc_joint *joint = &kc->segment[i - 1].joint;
        int joint_type = joint->type;

        // Inertia
        //

        // M_i^a = P_i^T M_i^A
        kcc_joint[joint_type].project_inertia(joint, &s.m_art[i], &s.m_app[i - 1]);

        // M_{i-1}^a' = {i-1}^X_i* M_i^a i^X_{i-1}
        mc_abi_tf_tgt_to_ref(&s.x_rel[i - 1], &s.m_app[i - 1], &s.m_tf[i]);

        // M_{i-1}^A += M_{i-1}^a'
        mc_abi_add(&s.m_art[i - 1], &s.m_tf[i], &s.m_art[i - 1]);


        // Force
        //

        // F_{bias,i}^A' = M_i^A Xdd_{bias,i}
        mc_abi_map_acc_twist_to_wrench(&s.m_art[i], &s.xdd_bias[i - 1], &s.f_bias_eom[i - 1]);

        // F_{bias,i}^A' += F_{bias,i}^A
        mc_wrench_add(&s.f_bias_eom[i - 1], &s.f_bias_art[i], &s.f_bias_eom[i - 1], 1);

        // F_{bias,i}^a = P_i^T F_{bias,i}^A'
        kcc_joint[joint_type].project_wrench(joint, &s.m_art[i], &s.f_bias_eom[i - 1], &s.f_bias_app[i - 1], 1);

        // F_{bias,i}^a' = {i-1}^X_i* F_{bias,i}^a
        mc_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_bias_app[i - 1], &s.f_bias_tf[i - 1], 1);

        // F_{bias,i-1}^A += F_{bias,i}^a'
        mc_wrench_add(&s.f_bias_art[i - 1], &s.f_bias_tf[i - 1], &s.f_bias_art[i - 1], 1);


        // F_{tau,i} = M_i^A S_i D^{-1} tau_{ff,i}
        kcc_joint[joint_type].ffd(joint, &s.m_art[i], &s.tau_ff[i - 1], &s.f_ff_jnt[i - 1], 1);

        // F_{tau,i}^a = P_i^T F_{tau,i}^A
        kcc_joint[joint_type].project_wrench(joint, &s.m_art[i], &s.f_ff_art[i], &s.f_ff_app[i - 1], 1);

        // F_{tau,i}^a += F_{tau,i}'
        mc_wrench_add(&s.f_ff_app[i - 1], &s.f_ff_jnt[i - 1], &s.f_ff_app[i - 1], 1);

        // F_{tau,i-1}^A = {i-1}^X_i* F_{tau,i}^a
        mc_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_ff_app[i - 1], &s.f_ff_art[i - 1], 1);


        // F_{ext,i}^a = P_i^T F_{ext,i}^A
        kcc_joint[joint_type].project_wrench(joint, &s.m_art[i], &s.f_ext_art[i], &s.f_ext_app[i - 1], 1);

        // F_{ext,i}^a' = {i-1}^X_i* F_{ext,i}^a
        mc_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_ext_app[i - 1], &s.f_ext_tf[i - 1], 1);

        // F_{ext,i-1}^A += F_{ext,i}^a'
        mc_wrench_add(&s.f_ext_art[i - 1], &s.f_ext_tf[i - 1], &s.f_ext_art[i - 1], 1);


        // F_{cstr,i}^a = P_i^T
        kcc_joint[joint_type].project_wrench(joint, &s.m_art[i], &s.f_cstr_art[i], &s.f_cstr_app[i - 1], NR_CSTR);

        // F_{cstr,i-1}^A = {i-1}^X_i* F_{cstr,i}^a
        mc_wrench_tf_tgt_to_ref(&s.x_rel[i - 1], &s.f_cstr_app[i - 1], &s.f_cstr_art[i - 1], NR_CSTR);


        // Energy
        //

        // E_{cstr,i-1}^A = E_{cstr,i}^A + (F_{cstr,i}^A)^T S_i D_i^{-1} S_i^T F_{cstr,i}^A
        kcc_joint[joint_type].decomp_e_cstr(joint,
                &s.m_art[i], &s.f_cstr_art[i],
                s.d_cstr_art[i], s.d_cstr_art[i - 1], NR_CSTR);
    }

    // nu_{cstr} = (E_{cstr,0})^{-1} E_{cstr,N}
    mc_eacc_balance(s.d_cstr_art[0], s.e_cstr, s.nu_cstr, NR_CSTR);

    mc_abi_log(&s.m_art[0]);
    mc_wrench_log(&s.f_bias_art[0], 1);
    mc_wrench_log(&s.f_ff_art[0], 1);
    mc_wrench_log(&s.f_ext_art[0], 1);
    mc_wrench_log(&s.f_cstr_art[0], NR_CSTR);


    for (int i = 1; i < s.nbody + 1; i++) {
        struct kcc_joint *joint = &kc->segment[i - 1].joint;
        int joint_type = joint->type;

        // Acceleration
        //

        // Xdd_{i-1}' = i^X_{i-1} Xdd_{i-1}
        gc_acc_twist_tf_ref_to_tgt(&s.x_rel[i - 1], &s.xdd[i - 1], &s.xdd_tf[i - 1]);

        // Xdd_{nact,i-1} = Xdd_{i-1}' + Xdd_{bias,i-1}
        gc_acc_twist_accumulate(&s.xdd_tf[i - 1], &s.xdd_bias[i - 1], &s.xdd_nact[i - 1]);


        // Force
        //

        // F_{nact,i} = M_i^A Xdd_{nact,i-1}'
        mc_abi_map_acc_twist_to_wrench(&s.m_art[i], &s.xdd_nact[i - 1], &s.f_bias_nact[i - 1]);

        // tau_{bias,i}^A = S^T F_{bias,i}
        kcc_joint[joint_type].ifk(joint, &s.f_bias_nact[i - 1], &s.tau_bias_art[i - 1], 1);


        // tau_{ff,i}^A = S^T F_{ff,i}
        kcc_joint[joint_type].ifk(joint, &s.f_ff_art[i], &s.tau_ff_art[i - 1], 1);


        // tau_{ext,i}^A = S^T F_{ext,i}
        kcc_joint[joint_type].ifk(joint, &s.f_ext_art[i], &s.tau_ext_art[i - 1], 1);


        // tau_{cstr,i}^A' = S^T F_{cstr,i}
        kcc_joint[joint_type].ifk(joint, &s.f_cstr_art[i], s.tau_cstr_art[i - 1], NR_CSTR);


        // Solve
        //
        int k = joint->revolute_joint.axis;
        double d = s.m_art[i].second_moment_of_mass.row[k].data[k] + joint->revolute_joint.inertia[0];

        s.tau_cstr[i - 1] = 0.0;
        for (int j = 0; j < NR_CSTR; j++) s.tau_cstr[i - 1] += s.nu_cstr[j] * s.tau_cstr_art[i - 1][j];

        s.tau_ctrl[i - 1] = s.tau_ff[i - 1] - s.tau_ff_art[i - 1] - s.tau_bias_art[i - 1] - s.tau_ext_art[i - 1] + s.tau_cstr[i - 1];
        s.qdd[i - 1] = s.tau_ctrl[i - 1] / d;


        // Resultant acceleration
        //

        // Xdd_{J,i} = S_i qdd_i
        kcc_joint[joint_type].fak(joint, &s.qdd[i - 1], &s.xdd_jnt[i - 1]);

        // Xdd_{J,i}' = Xdd_{J,i} + Xdd_{bias,i}
        gc_acc_twist_add(&s.xdd_jnt[i - 1], &s.xdd_bias[i - 1], &s.xdd_net[i - 1]);

        // Xdd_i = Xdd_{i-1}' + Xdd_{J,i}'
        gc_acc_twist_accumulate(&s.xdd_tf[i - 1], &s.xdd_net[i - 1], &s.xdd[i]);
    }

    gc_acc_twist_log(&s.xdd[s.nbody]);
}


int main(int argc, char **argv)
{
    achd_a();
    achd_c();

    return 0;
}
