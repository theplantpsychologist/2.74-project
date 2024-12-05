function output = simulate_twolegs(slope, angular_velocity)
    %% Define fixed paramters
    m1 =.0393 + .2;         m2 =.0368; 
    m3 = .00783;            m4 = .0155;
    I1 = 25.1 * 10^-6;      I2 = 53.5 * 10^-6;
    I3 = 9.25 * 10^-6;      I4 = 22.176 * 10^-6;
    l_OA=.011;              l_OB=.042; 
    l_AC=.096;              l_DE=.091;
    l_O_m1=0.032;           l_B_m2=0.0344; 
    l_A_m3=0.0622;          l_C_m4=0.0610;
    N = 18.75;
    Ir = 0.0035/N^2;
    g = 9.81;    

    q1_min = deg2rad(-90); %-80
    q1_max = deg2rad(90); %80
    q2_min = deg2rad(25); %10
    q2_max = deg2rad(135); %170
    qlims = [q1_min, q1_max, q2_min, q2_max];

    %motor parameters
    Kt = 0.191; % torque constant of motor
    R_internal = 1.594; % internal resistance (mbed lab 3) 
    v = 12; %volts

    restitution_coeff = 0.;
    friction_coeff = 0.8;
    ground_height = 0;

    step_th = slope;
    % step_w = 0.05; % step depth
    % step_h = tan(step_th)*step_w; %find step height. tan(stair_th) = stair_h/stair_w

    step_h = 1*0.0254;
    step_w = step_h/tan(step_th);

    % offset_x = -step_w/3;
    offset_x = 0;

    %%%%%%%%%%%%%%%%%%% CHANGE GROUND CONSTRAINT HERE %%%%%%%%%%%%%%%%%%
    stair_height = @(x) floor((x+offset_x)./step_w).*step_h + ground_height;

    %% Parameter vector
    p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g R_internal Kt v]';
       
    % %% Simulation Parameters Set 2 -- Operational Space Control

    %%%%%%%%%%%%%% CONTROL ELLIPSE PARAMETERS %%%%%%%%%%%%%%
    p_traj.omega = -angular_velocity; %rad/sec
    L = (l_DE+l_OB);
    p_traj.x_off   = -L*0.5;
    p_traj.y_off   = -L*1.3;
    p_traj.a     = 0.08;
    p_traj.b     = 0.06;
    p_traj.phi  = slope+ deg2rad(20);
    p_traj.phase = pi;


    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 5;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    % Initial state: [x; y; th1; th2; th3; th4; dx; dy; dth1; dth2; dth3; dth4]
    z0 = [0.0692; 0.1868; -pi/4; pi/3; -pi/4; pi/3; 0; 0; 0; 0; 0; 0];
    z_out = zeros(12,num_step);
    z_out(:,1) = z0;
    
    torque_comm = zeros(6,num_step); 
    actual_torques = zeros(4, num_step);

    for i=1:num_step-1
        
        dz = dynamics(tspan(i), z_out(:,i), p, p_traj);
        % Update state with dynamics
        z_temp = z_out(:,i) + dz * dt;
    
        % Check for ground contact and update velocities if necessary
        z_temp(7:12) = discrete_impact_contact(z_temp, p,restitution_coeff,friction_coeff,stair_height);
    
        % Check for joint limit constraint and update velocities if necessary
        z_temp(7:12) = joint_limit_constraint(z_temp, p, qlims);
    
        % Update the state
        z_out(:,i+1) = z_temp;

        % Updates new column of torque_comm at each time step 
        torque_comm(:,i) = control_law(tspan(i), z_out(:,i), p, p_traj);

        % get clipped torques
        actual_torques(:,i) = motorDynamics(tspan(i),z_out(:,i),p,torque_comm(:,i));
       
    end
    totalEnergy = computeMotorEnergySpent(tspan, torque_comm, Kt, R_internal);
    output = [z_out(2,end), totalEnergy(end), slope, angular_velocity];
    %% Compute Energy
    % E = potential_energy_twolegs(z_out,p);
    % figure(1); clf
    % plot(tspan,E); xlabel('Time (s)'); ylabel('Energy (J)');
    % 
    % %% Compute foot position over time
    % rE1 = zeros(2,length(tspan));
    % vE1 = zeros(2,length(tspan));
    % rE2 = zeros(2,length(tspan));
    % vE2 = zeros(2,length(tspan));
    % for i = 1:length(tspan)
    %     rE1(:,i) = position_foot1(z_out(:,i),p);
    %     vE1(:,i) = velocity_foot1(z_out(:,i),p);
    %     rE2(:,i) = position_foot2(z_out(:,i),p);
    %     vE2(:,i) = velocity_foot2(z_out(:,i),p);
    % end
    % 
    % figure(2); clf;
    % plot(tspan,rE1(1,:),'r','LineWidth',2)
    % hold on
    % % plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    % plot(tspan,rE1(2,:),'b','LineWidth',2)
    % % plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
    % plot(tspan,rE2(1,:),'m','LineWidth',2)
    % plot(tspan,rE2(2,:),'c','LineWidth',2)
    % hold off %added 
    % 
    % xlabel('Time (s)'); ylabel('Position (m)'); legend({'x1','x2','y1','y2'});
    % % xlabel('Time (s)'); ylabel('Position (m)'); legend({'x1','x_d','y1','y_d'});
    % figure(3); clf;
    % plot(tspan,vE1(1,:),'r','LineWidth',2)
    % hold on
    % plot(tspan,vE1(2,:),'b','LineWidth',2)
    % hold on
    % plot(tspan,vE2(1,:),'m','LineWidth',2)
    % hold on
    % plot(tspan,vE2(2,:),'c','LineWidth',2)
    % 
    % xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x1','vel_y1','vel_x2','vel_y2'});
    % 
    % figure(4)
    % plot(tspan,z_out(1:4,:)*180/pi) %convert gen coord: rad to degrees
    % legend('q1','q2','q3','q4');
    % xlabel('Time (s)');
    % ylabel('Angle (deg)');
    % 
    % figure(5)
    % plot(tspan,z_out(5:8,:)*180/pi) %convert gen vel's: rad/s to degrees/s
    % legend('q1dot','q2dot','q3dot','q4dot');
    % xlabel('Time (s)');
    % ylabel('Angular Velocity (deg/sec)');
    % 
    % %% Animate Solution
    % figure(6); clf;
    % hold on
    % 
    % % Ground Q2.3
    % steps = 50;
    % stair_X = linspace(-step_w,step_w*steps,1000);
    % stair_Y = stair_height(stair_X);
    % plot(stair_X,stair_Y,'k');
    % % plot([-2 5], [ground_height ground_height], 'k--');
    % 
    % animateSol(tspan, z_out,p, p_traj);
    % 
    % % Calculate & Plot Energy spent by motors (based on currents and torques commanded)
    % figure(7); clf 
    % totalEnergy = computeMotorEnergySpent(tspan, torque_comm, Kt, R_internal);
    % net_Energy = totalEnergy - E;
    % % motorEnergy = tspan;
    % plot(tspan,net_Energy); hold on;
    % plot(tspan,E);
    % plot(tspan,totalEnergy);
    % legend('dissipated','potential','motor');
    % xlabel('Time (s)'); ylabel('Energy (J)');
    % 
    % figure(8); clf
    % plot(tspan, torque_comm(3,:)); hold on; 
    % % plot(tspan, torque_comm(4,:)); hold on;
    % % plot(tspan, torque_comm(5,:)); hold on;
    % % plot(tspan, torque_comm(6,:)); hold on;
    % plot(tspan, actual_torques(1,:));
    % % plot(tspan, actual_torques(2,:));
    % % plot(tspan, actual_torques(3,:));
    % % plot(tspan, actual_torques(4,:));
    % legend('LEG1: Motor 1','LEG1: Actual 1')
    % xlabel('Time (s)'); ylabel('Torque Output by Motor (N*m)');
    % 
    % figure(9); clf
    % plot(tspan,torque_comm(3,:)./Kt)
    % xlabel('Time (s)'); ylabel('Current(Amp)');

end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 700.; % Spring stiffness X
    K_y = 700.; % Spring stiffness Y
    D_x = 10.;  % Damping X
    D_y = 10.;  % Damping Y

    K = [K_x, 0 ; 0, K_y]; %K_xy = 0 
    D = [D_x, 0 ; 0, D_y]; %D_xy = 0 
    
    %base position
    x_center = z(1) + p_traj.x_off;
    y_center = z(2) + p_traj.y_off;
    % 
    %current ellipse angles
    phase = p_traj.phase; %phase shift between legs
    w = p_traj.omega;
    a = p_traj.a;
    b = p_traj.b;
    phi = p_traj.phi;
    theta1 = w*t; %foot 1
    theta2 = theta1 + phase;

    %ellipse equation
    x_pos = @(theta) x_center + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi);
    y_pos = @(theta) y_center + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi);

    x_vel = @(theta) - a * sin(theta) * cos(phi)*w - b * cos(theta) * sin(phi) * w;
    y_vel = @(theta) - a * sin(theta) * sin(phi)*w + b * cos(theta) * cos(phi) * w;

    x_acc = @(theta) - a * cos(theta) * cos(phi)*w^2 + b * sin(theta) * sin(phi) * w^2;
    y_acc = @(theta) - a * cos(theta) * sin(phi)*w^2 - b * sin(theta) * cos(phi) * w^2;

    % Desired position of foot is a circle
    rEd1 = [x_pos(theta1) y_pos(theta1)]';
    vEd1 = [x_vel(theta1) y_vel(theta1)]';
    aEd1 = [x_acc(theta1) y_acc(theta1)]';
    % Leg 2 
    rEd2 = [x_pos(theta2) y_pos(theta2)]';
    vEd2 = [x_vel(theta2) y_vel(theta2)]';
    aEd2 = [x_acc(theta2) y_acc(theta2)]';

    % Actual position and velocity 
    rE1 = position_foot1(z,p);
    vE1 = velocity_foot1(z,p);
    rE2 = position_foot2(z,p);
    vE2 = velocity_foot2(z,p);


    A  = A_twolegs(z,p);       % 4x4 matrix,now 6x6?
    inv_A = inv(A);
    Corr = Corr_twolegs(z,p);      % 6x2 matrix
    Grav = Grav_twolegs(z,p);      % 6x2 matrix
    qdot = qdot_twolegs(z,p);  % 6x1 matrix
    J1  = jacobian_foot1(z,p); % 2x6 matrix
    J2  = jacobian_foot2(z,p); % 2x6 matrix

    dJ1 = jacobian_dot_foot1(z,p);
    dJ2 = jacobian_dot_foot2(z,p);
    % Mass matrix of actuated joints
    % Compute Jacobians with respect to actuated joints only
    J1 = jacobian_foot1(z, p); % Original 2x6 Jacobian
    J2 = jacobian_foot2(z, p); % Original 2x6 Jacobian

    % To get f, WE NEED LAMBDA, MU, & RHO
    lambda = inv (J1 * inv_A * (J1.'));
    mu1 = (lambda * J1 * inv_A * Corr) - (lambda * dJ1 * qdot);
    rho1 = lambda * J1 * inv_A * Grav;

    % Compute virtual foce 
    intermediate1 = aEd1 + K * (rEd1 - rE1(1:2, 1:end)) + D * (vEd1 - vE1(1:2, 1:end));
    f1 = lambda * intermediate1 + mu1 + rho1; % should result in a 2x1 matrix

    lambda2 = inv (J2 * inv_A * (J2.'));
    mu2 = (lambda2 * J2 * inv_A * Corr) - (lambda2 * dJ2 * qdot);
    rho2 = lambda2 * J2 * inv_A * Grav;
    intermediate2 = aEd2+ K * (rEd2 - rE2(1:2, 1:end)) + D * (vEd2 - vE2(1:2, 1:end)); % aEd2
    f2 = lambda2 * intermediate2 + mu2 + rho2;


    %% Task-space compensation and feed forward for Question 1.8

    % Map to joint torques  

    tau1 = J1' * f1;
    tau1 = tau1(3:4); 

    tau2 = J2' * f2;
    tau2 = tau2(5:6);
    tau = [0;0;tau1; tau2]; %% is a 6x1 matrix
end

function dz = dynamics(t,z,p, p_traj)
    % Get mass matrix
    A = A_twolegs(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);

    % Get b = Q - V(q,qd) - G(q)
    b = b_twolegs(z,tau,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = zeros(size(z));
    
    dz(1:6) = z(7:12);
    dz(7:12) = qdd;
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC_fun)

    % Actual position and velocity 
    rE1 = position_foot1(z,p);
    vE1 = velocity_foot1(z,p);
    rE2 = position_foot2(z,p);
    vE2 = velocity_foot2(z,p);
    J1  = jacobian_foot1(z,p); % 2x6 matrix
    J2  = jacobian_foot2(z,p); % 2x6 matrix
    A  = A_twolegs(z,p); %6x6 matrix
    inv_A = inv(A);
    osim1 = inv (J1 * inv_A * J1');
    osim2 = inv (J2 * inv_A * J2');
    xhat = [1,0];
    yhat = [0,1];
    mass_eff_x1 = 1 / (xhat * inv(osim1) * xhat'); %m_eff_r = 1/(rhat*inv(osim)*rhat.')
    mass_eff_y1 = 1 / (yhat * inv(osim1) * yhat');
    mass_eff_x2 = 1 / (xhat * inv(osim2) * xhat'); %m_eff_r = 1/(rhat*inv(osim)*rhat.')
    mass_eff_y2 = 1 / (yhat * inv(osim2) * yhat');

    x1 = rE1(1);
    x2 = rE2(1);
    yC1 = yC_fun(x1);
    yC2 = yC_fun(x2);

    C_y1 = rE1(2) - yC1;
    Cdot_y1 = vE1(2);

    C_y2 = rE2(2) - yC2;
    Cdot_y2 = vE2(2);

    qdot = qdot_twolegs(z,p); % 6x1, to be updated
    % Leg 1
    if (C_y1 <= 0) && (Cdot_y1 <= 0) 
        F_impulse_y1 = mass_eff_y1 * (-rest_coeff * Cdot_y1 - J1(2,:) * qdot);
        F_impulse_x1 = mass_eff_x1 * (-J1(1,:) * qdot);
        if abs(F_impulse_x1) > abs(fric_coeff * F_impulse_y1)
            F_impulse_x1 = sign(F_impulse_x1) * fric_coeff * F_impulse_y1;
        end
        %update both for x and y impulses
        qdot = qdot + (inv_A * J1(2,:).' * F_impulse_y1) + (inv_A * J1(1,:).' * F_impulse_x1); %update qdot1 & qdot2
    end
    % Leg 2
    if (C_y2 <= 0) && (Cdot_y2 <= 0) 
        F_impulse_y2 = mass_eff_y2 * (-rest_coeff * Cdot_y2 - J2(2,:) * qdot);
        F_impulse_x2 = mass_eff_x2 * (-J2(1,:) * qdot);
        if abs(F_impulse_x2) > abs(fric_coeff * F_impulse_y2)
            F_impulse_x2 = sign(F_impulse_x2) * fric_coeff * F_impulse_y2;
        end
        %update both for x and y impulses
        qdot = qdot + (inv_A * J2(2,:).' * F_impulse_y2) + (inv_A * J2(1,:).' * F_impulse_x2); %update qdot3 & qdot4
    end
    qdot = qdot;

end

function totalEnergy = computeMotorEnergySpent (t, torque_comm, Kt, R_internal)
    % calculates energy consumed by motors
    % returns totalEnergy vector (useful for plotting) 

    torque_comm = torque_comm(3:6, :); % entries 1&2 outputted by control_law aren't motor torques 
    totalEnergy = zeros(1,length(t));
    % totalEnergy = 0;
    power = zeros(1,length(t));
    % Calculate the enrgy motor by motor, then add to totalEnergy
    for m = 1:4

        % Calculate power vector & integrate at each step to find energy up
        % to your current stage, then update totalEnergy at that stage
        for c = 1:length(t)
            current = torque_comm(m,c)/Kt;
            power(c) = (current.^2) / R_internal;
            % trapz over the entire time and power vectors works because
            % the power vector is zero beyond the current stage 
            energyTemp = trapz(t, power);
            totalEnergy(c) = totalEnergy(c) + energyTemp;
        end
        % Reset power vector for the power & energy calculations of the next motor (using next row of torque_comm)
        power = zeros(1,length(t)); 
    end

    totalEnergy = totalEnergy;
end

function actual_torques = motorDynamics (t,z,p,torque_comm)
    tau_comm = torque_comm(3:6); 
    v_max = p(22); Kt = p(21); R = p(20); %motor parameters
    w = z(9:12); % [w1, w2, w3, w4]
    I_max = v_max/R; %assumes no back-emf
    stall_torque = abs(I_max*Kt);
    back_emf = abs(Kt^2/R.*w);
    Max_torques = stall_torque - back_emf; %find max torque speed curve at v_max
    actual_torques = min(Max_torques, abs(tau_comm)).*sign(tau_comm);%apply saturation limit
end

function qdot = joint_limit_constraint(z, p, limits)
    % Extract positions and velocities and limits
    q = z(1:6);      % Positions: [x; y; th1; th2; th3; th4]
    dq = z(7:12);     % Velocities: [dx; dy; dth1; dth2; dth3; dth4]
    qdot = dq;       % Initialize qdot with current velocities
    q1min = limits(1);
    q1max = limits(2);
    q2min = limits(3);
    q2max = limits(4);

    q1 = q(3);       % Joint angle th1
    q2 = q(4);       % Joint angle th2
    q3 = q(5);
    q4 = q(6);
    dth1 = dq(3);    % Joint velocity dth1
    dth2 = dq(4);    % Joint velocity dth2
    dth3 = dq(5);
    dth4 = dq(6);
    
    A = A_twolegs(z, p);% Mass matrix
    inv_mass = inv(A);

    if ((q1 <= q1min) && (dth1 < 0)) || ((q1 >= q1max) && (dth1 > 0))% Joint limit violated and moving further into violation
        Jc = [0, 0, 1, 0, 0, 0]; % Compute the Jacobian of the constraint c(q) = q1 - q1_min >= 0
        
        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc'); % -(1 + e) * (Jc * dq) is equivalent to -(gamma * dCy + Jc * dq) because ddtCy = dCy/q*dq = Jc*dq
        qdot = qdot + inv_mass * (Jc' * Lambda);
    end

    if ((q2 <= q2min) && (dth2 < 0)) || ((q2 >= q2max) && (dth2 > 0))
        Jc = [0, 0, 0, 1, 0, 0];% c(q) = q2 - q2_max <= 0
        A = A_twolegs(z, p); % Mass matrix
        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc');% impulse force
        qdot = qdot + inv_mass * (Jc' * Lambda); % Adjust velocities
    end
    if ((q3 <= q1min) && (dth3 < 0)) || ((q3 >= q1max) && (dth3 > 0))
        Jc = [0, 0, 0, 0, 1, 0];% c(q) = q2 - q2_max <= 0
        A = A_twolegs(z, p); % Mass matrix
        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc');% impulse force
        qdot = qdot + inv_mass * (Jc' * Lambda); % Adjust velocities
    end
    if ((q4 <= q2min) && (dth4 < 0)) || ((q4 >= q2max) && (dth4 > 0))
        Jc = [0, 0, 0, 0, 0, 1];% c(q) = q2 - q2_max <= 0
        A = A_twolegs(z, p); % Mass matrix
        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc');% impulse force
        qdot = qdot + inv_mass * (Jc' * Lambda); % Adjust velocities
    end
end


function animateSol(tspan,x,p,p_traj)
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
    h_OB2 = plot([0],[0],'LineWidth',2);
    h_AC2 = plot([0],[0],'LineWidth',2);
    h_BD2 = plot([0],[0],'LineWidth',2);
    h_CE2 = plot([0],[0],'LineWidth',2);
    h_ellipse = plot([0],[0],'LineWidth',2);
    h_ref1 = plot([0],[0],'r.');
    h_ref2 = plot([0],[0],'b.');
   
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-1 2.0 -0.25 1.5]);

    TH = 0:.1:2*pi;
    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_twolegs(z,p);
        
        rO = keypoints(:,1); % Vector to base of cart
        rA = keypoints(:,2);
        rB = keypoints(:,3); % Vector to tip of pendulum
        rC = keypoints(:,4);
        rD = keypoints(:,5);
        rE = keypoints(:,6); % Vector to base of cart
        rA2 = keypoints(:,7);
        rB2 = keypoints(:,8); % Vector to tip of pendulum
        rC2 = keypoints(:,9);
        rD2 = keypoints(:,10);
        rE2 = keypoints(:,11);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[rO(1) rB(1)]);
        set(h_OB,'YData',[rO(2) rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);

        set(h_OB2,'XData',[rO(1) rB2(1)]);
        set(h_OB2,'YData',[rO(2) rB2(2)]);
        
        set(h_AC2,'XData',[rA2(1) rC2(1)]);
        set(h_AC2,'YData',[rA2(2) rC2(2)]);
        
        set(h_BD2,'XData',[rB2(1) rD2(1)]);
        set(h_BD2,'YData',[rB2(2) rD2(2)]);
        
        set(h_CE2,'XData',[rC2(1) rE2(1)]);
        set(h_CE2,'YData',[rC2(2) rE2(2)]);

        %draw ellipse
        x_center = z(1) + p_traj.x_off;
        y_center = z(2) + p_traj.y_off;
        X_ell = @(theta) x_center + p_traj.a .* cos(theta) .* cos(p_traj.phi) - p_traj.b .* sin(theta) .* sin(p_traj.phi);
        Y_ell = @(theta) y_center + p_traj.a * cos(theta) .* sin(p_traj.phi) + p_traj.b .* sin(theta) .* cos(p_traj.phi);

        set(h_ellipse, 'XData',X_ell(TH));
        set(h_ellipse, 'YData',Y_ell(TH));
            
        %track references
        theta1 = t*p_traj.omega;
        theta2 = t*p_traj.omega+p_traj.phase;
        set(h_ref1, {'XData','YData'}, {X_ell(theta1), Y_ell(theta1)});
        set(h_ref2, {'XData','YData'}, {X_ell(theta2), Y_ell(theta2)});
        

        pause(.01)
    end
end