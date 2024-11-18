function [] = buildTwoLegDynamics()
    % q         = [x; y; th1  ; th2 ; th3 ; th4 ];      % generalized coordinates
    % q_dot     = [dx; dy; dth1 ; dth2 ; dth3 ; dth4];    % first time derivatives
    % q_ddot    = [ddx; ddy; ddth1;ddth2;ddth3;ddth4];  % second time derivatives
    % u         = [Fx; Fy; tau1 ; tau2; tau3 ; tau4];     % controls
    % not used: F         = [Fx ; Fy];

    % Parameters
    % p         = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';

    q = casadi.SX.sym('q', 6, 1); % I'm guessing the # after the var is the size 
    q_dot = casadi.SX.sym('q_dot', 6, 1);
    q_ddot = casadi.SX.sym('q_ddot', 6, 1);
    u = casadi.SX.sym('u', 6, 1);
    %F = casadi.SX.sym('F', 2, 1);
    par = casadi.SX.sym('p', 19, 1);
    par_mass = par(1:4);
    par_inertia = par(5:8);
    par_inertia_rotor = par(9);
    par_gear_ratio = par(10);
    par_com_distances = par(11:14);
    par_lengths = par(15:18);
    par_gravity = par(19);

    % i_hat = [0; -1; 0]; % old version 
    % j_hat = [1; 0; 0];
    % e1_hat = cos(q(1))*i_hat + sin(q(1))*j_hat;
    % e2_hat = cos(q(1) + q(2))*i_hat + sin(q(1) + q(2))*j_hat;
    ihat = [1; 0; 0]; %same as x_hat
    jhat = [0; 1; 0]; %same as y_hat
    k_hat = cross(i_hat, j_hat);
  
    ehat1 =  sin(q(3))*ihat - cos(q(3))*jhat;
    ehat2 =  sin(q(3)+q(4))*ihat - cos(q(3)+q(4))*jhat;
    ehat3 =  sin(q(5))*ihat - cos(q(5))*jhat;
    ehat4 =  sin(q(5)+q(6))*ihat - cos(q(5)+q(6))*jhat;

    rO = q(1)*ihat + q(2)*jhat;
    rOA1 = rO + par_lengths(1) * ehat1;  
    rOB1 = rO + par_lengths(2) * ehat1;
    rOC1 = rOA1 + par_lengths(3) * ehat2; 
    rOD1 = rOB1 + par_lengths(3) * ehat2; 
    rOE1 = rOD1 + par_lengths(4) * ehat1; 
    rOA2 = rO + par_lengths(1) * ehat3;  
    rOB2 = rO + par_lengths(2) * ehat3;
    rOC2 = rOA2 + par_lengths(3) * ehat4; 
    rOD2 = rOB2 + par_lengths(3) * ehat4; 
    rOE2 = rOD2 + par_lengths(4) * ehat3; 
    % r_A = par_lengths(1) * e1_hat;
    % r_B = par_lengths(2) * e1_hat;
    % r_C = r_A + par_lengths(3) * e2_hat;
    % r_D = r_B + par_lengths(3) * e2_hat;
    % r_E = r_D + par_lengths(4) * e1_hat;

    r_m1_1 = rO + par_com_distances(1) * ehat1; 
    r_m2_1 = rOB1 + par_com_distances(2) * ehat2;
    r_m3_1 = rOA1 + par_com_distances(3) * ehat2;
    r_m4_1 = rOC1 + par_com_distances(4) * ehat1;
    r_m1_2 = rO + par_com_distances(1) * ehat3; 
    r_m2_2 = rOB2 + par_com_distances(2) * ehat4;
    r_m3_2 = rOA2 + par_com_distances(3) * ehat4;
    r_m4_2 = rOC2 + par_com_distances(4) * ehat3;
    % r_m1 = par_com_distances(1) * e1_hat;
    % r_m2 = r_B + par_com_distances(2) * e2_hat;
    % r_m3 = r_A + par_com_distances(3) * e2_hat;
    % r_m4 = r_C + par_com_distances(4) * e1_hat;

    r_A_dot1 = casadi_ddt(rOA1, [q; q_dot], [q_dot; q_ddot]);
    r_B_dot1 = casadi_ddt(rOB1, [q; q_dot], [q_dot; q_ddot]);
    r_C_dot1 = casadi_ddt(rOC1, [q; q_dot], [q_dot; q_ddot]);
    r_D_dot1 = casadi_ddt(rOD1, [q; q_dot], [q_dot; q_ddot]);
    r_E_dot1 = casadi_ddt(rOE1, [q; q_dot], [q_dot; q_ddot]);
    r_A_dot2 = casadi_ddt(rOA2, [q; q_dot], [q_dot; q_ddot]);
    r_B_dot2 = casadi_ddt(rOB2, [q; q_dot], [q_dot; q_ddot]);
    r_C_dot2 = casadi_ddt(rOC2, [q; q_dot], [q_dot; q_ddot]);
    r_D_dot2 = casadi_ddt(rOD2, [q; q_dot], [q_dot; q_ddot]);
    r_E_dot2 = casadi_ddt(rOE2, [q; q_dot], [q_dot; q_ddot]);
    % r_A_dot = casadi_ddt(r_A, [q; q_dot], [q_dot; q_ddot]);
    % r_B_dot = casadi_ddt(r_B, [q; q_dot], [q_dot; q_ddot]);
    % r_C_dot = casadi_ddt(r_C, [q; q_dot], [q_dot; q_ddot]);
    % r_D_dot = casadi_ddt(r_D, [q; q_dot], [q_dot; q_ddot]);
    % r_E_dot = casadi_ddt(r_E, [q; q_dot], [q_dot; q_ddot]);

    r_m1_dot1 = casadi_ddt(r_m1_1, [q; q_dot], [q_dot; q_ddot]);
    r_m2_dot1 = casadi_ddt(r_m2_1, [q; q_dot], [q_dot; q_ddot]);
    r_m3_dot1 = casadi_ddt(r_m3_1, [q; q_dot], [q_dot; q_ddot]);
    r_m4_dot1 = casadi_ddt(r_m4_1, [q; q_dot], [q_dot; q_ddot]);
    r_m1_dot2 = casadi_ddt(r_m1_2, [q; q_dot], [q_dot; q_ddot]);
    r_m2_dot2 = casadi_ddt(r_m2_2, [q; q_dot], [q_dot; q_ddot]);
    r_m3_dot2 = casadi_ddt(r_m3_2, [q; q_dot], [q_dot; q_ddot]);
    r_m4_dot2 = casadi_ddt(r_m4_2, [q; q_dot], [q_dot; q_ddot]);
    % r_m1_dot = casadi_ddt(r_m1, [q; q_dot], [q_dot; q_ddot]);
    % r_m2_dot = casadi_ddt(r_m2, [q; q_dot], [q_dot; q_ddot]);
    % r_m3_dot = casadi_ddt(r_m3, [q; q_dot], [q_dot; q_ddot]);
    % r_m4_dot = casadi_ddt(r_m4, [q; q_dot], [q_dot; q_ddot]);

    w_1 = q_dot(3); % dth1
    w_1and2 = q_dot(3) + q_dot(4); % dth1 + dth2
    w_3 = q_dot(5); % dth3
    w_3and4 = q_dot(5) + q_dot(6); % dth3 + dth4

    % Kinetic energy of masses in motion: translational + rotational
    T1_1 = 0.5 * par_mass(1) * dot(r_m1_dot1, r_m1_dot1) + 0.5 * par_inertia(1) * w_1^2;
    T1_2 = 0.5 * par_mass(2) * dot(r_m2_dot1, r_m2_dot1) + 0.5 * par_inertia(2) * w_1and2^2;
    T1_3 = 0.5 * par_mass(3) * dot(r_m3_dot1, r_m3_dot1) + 0.5 * par_inertia(3) * w_1and2^2;
    T1_4 = 0.5 * par_mass(4) * dot(r_m4_dot1, r_m4_dot1) + 0.5 * par_inertia(4) * w_1^2;
    T2_1 = 0.5 * par_mass(1) * dot(r_m1_dot2, r_m1_dot2) + 0.5 * par_inertia(1) * w_3^2;
    T2_2 = 0.5 * par_mass(2) * dot(r_m2_dot2, r_m2_dot2) + 0.5 * par_inertia(2) * w_3and4^2;
    T2_3 = 0.5 * par_mass(3) * dot(r_m3_dot2, r_m3_dot2) + 0.5 * par_inertia(3) * w_3and4^2;
    T2_4 = 0.5 * par_mass(4) * dot(r_m4_dot2, r_m4_dot2) + 0.5 * par_inertia(4) * w_3^2;
    % Gears
    T1_1_rotor = 0.5 * par_inertia_rotor * (par_gear_ratio * q_dot(3))^2;
    T1_2_rotor = 0.5 * par_inertia_rotor * (q_dot(3) + par_gear_ratio * q_dot(4))^2;
    T2_1_rotor = 0.5 * par_inertia_rotor * (par_gear_ratio * q_dot(5))^2;
    T2_2_rotor = 0.5 * par_inertia_rotor * (q_dot(5) + par_gear_ratio * q_dot(6))^2;

    V1_g1 = par_mass(1) * par_gravity * dot(r_m1_1, j_hat);
    V1_g2 = par_mass(2) * par_gravity * dot(r_m2_1, j_hat);
    V1_g3 = par_mass(3) * par_gravity * dot(r_m3_1, j_hat);
    V1_g4 = par_mass(4) * par_gravity * dot(r_m4_1, j_hat);
    V2_g1 = par_mass(1) * par_gravity * dot(r_m1_2, j_hat);
    V2_g2 = par_mass(2) * par_gravity * dot(r_m2_2, j_hat);
    V2_g3 = par_mass(3) * par_gravity * dot(r_m3_2, j_hat);
    V2_g4 = par_mass(4) * par_gravity * dot(r_m4_2, j_hat);

    T = simplify(T1_1 + T1_2 + T1_3 + T1_4 + T2_1 + T2_2 + T2_3 + T2_4 + T1_1_rotor + T1_2_rotor + T2_1_rotor + T2_2_rotor);
    V = simplify(V1_g1 + V1_g2 + V1_g3 + V1_g4 + V2_g1 + V2_g2 + V2_g3 + V2_g4);

    Q_tau1 = casadi_M2Q(u(3) * k_hat, w_1 * k_hat, q_dot);
    Q_tau2 = casadi_M2Q(u(4) * k_hat, w_2 * k_hat, q_dot);
    Q_tau2_reaction = casadi_M2Q(-u(4) * k_hat, w_1 * k_hat, q_dot);
    Q_tau3 = casadi_M2Q(u(5) * k_hat, w_3 * k_hat, q_dot);
    Q_tau4 = casadi_M2Q(u(6) * k_hat, w_3and4 * k_hat, q_dot);
    Q_tau4_reaction = casadi_M2Q(-u(6) * k_hat, w_3 * k_hat, q_dot);
    Q_tau = simplify(Q_tau1 + Q_tau2 + Q_tau3 + Q_tau4 + Q_tau2_reaction + Q_tau4_reaction);

    Q_force_x = casadi_F2Q(F(1)*ihat, rO * ihat);
    Q_force_y = casadi_F2Q(F(2)*jhat, rO * jhat);
    Q_force = Q_force_x + Q_force_y;
    
    Q = Q_tau + Q_force;

    keypoints = [rO(1:2) rOA1(1:2) rOB1(1:2) rOC1(1:2) rOD1(1:2) rOE1(1:2) rOA2(1:2) rOB2(1:2) rOC2(1:2) rOD2(1:2) rOE2(1:2)];

    E = simplify(T + V);
    L = simplify(T - V);
    EoM = casadi_ddt(jacobian(L, q_dot).', [q; q_dot], [q_dot; q_ddot]) - jacobian(L, q).' - Q;
    EoM_leg = casadi.Function('EoM_leg', {q, q_dot, q_ddot, u, par}, {EoM});

    A = (jacobian(EoM, q_ddot));
    b = -EoM_leg(q, q_dot, 0, u, par); % unsure if I have to pass F as well
    J1 = jacobian(rOE1, q); J1 = J1(1:2, 1:6);
    J2 = jacobian(rOE2, q); J2 = J2(1:2, 1:6);
    z = [q; q_dot];

    A_fn = casadi.Function('A_fn', {z, par}, {A});
    b_fn = casadi.Function('b_fn', {z, u, par}, {b});
    energy_fn = casadi.Function('energy_fn', {z, par}, {E});
    pos_end_effector1 = casadi.Function('pos_end_effector1', {z, par}, {rOE1(1:2)}); 
    vel_end_effector1 = casadi.Function('vel_end_effector1', {z, par}, {r_E_dot1(1:2)});
    J_end_effector1 = casadi.Function('J_end_effector1', {z, par}, {J1});
    pos_end_effector2 = casadi.Function('pos_end_effector2', {z, par}, {rOE2(1:2)});
    vel_end_effector2 = casadi.Function('vel_end_effector2', {z, par}, {r_E_dot2(1:2)});
    J_end_effector2 = casadi.Function('J_end_effector2', {z, par}, {J2});
    keypoints_fn = casadi.Function('keypoints_fn', {z, par}, {keypoints});
    
    A_fn.save('codegen/A_fn.casadi');
    b_fn.save('codegen/b_fn.casadi');
    energy_fn.save('codegen/energy_fn.casadi');
    pos_end_effector1.save('codegen/pos_end_effector1.casadi');
    vel_end_effector1.save('codegen/vel_end_effector1.casadi');
    J_end_effector1.save('codegen/J_end_effector1.casadi');
    pos_end_effector2.save('codegen/pos_end_effector2.casadi');
    vel_end_effector2.save('codegen/vel_end_effector2.casadi');
    J_end_effector2.save('codegen/J_end_effector2.casadi');
    keypoints_fn.save('codegen/keypoints_fn.casadi');
end