function simulate_leg()
    %% Definte fixed paramters
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
    
    restitution_coeff = 0.;
    friction_coeff = 0.3;
    ground_height = -0.13;
    %% Parameter vector
    p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';
       
    %% Simulation Parameters Set 2 -- Operational Space Control
    p_traj.omega = 3; %rad/sec
    p_traj.x_0   = 0;
    p_traj.y_0   = -.125;
    p_traj.r     = 0.025;
    
    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 5;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 
    z0 = [-pi/4; pi/2; 0; 0];
    z_out = zeros(4,num_step);
    z_out(:,1) = z0;
    
    for i=1:num_step-1
        dz = dynamics(tspan(i), z_out(:,i), p, p_traj);
        % Velocity update with dynamics
        z_out(:,i+1) = z_out(:,i) + dz*dt;
        
        % Position update
        z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt;
    end
    
    %% Compute Energy
    E = energy_leg(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    
    %% Compute foot position over time
    rE1 = zeros(2,length(tspan));
    vE1 = zeros(2,length(tspan));
    rE2 = zeros(2,length(tspan));
    vE2 = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rE1(:,i) = position_foot1(z_out(:,i),p);
        vE1(:,i) = velocity_foot1(z_out(:,i),p);
        rE2(:,i) = position_foot2(z_out(:,i),p);
        vE2(:,i) = velocity_foot2(z_out(:,i),p);
    end
    
    figure(2); clf;
    plot(tspan,rE1(1,:),'r','LineWidth',2)
    hold on
    % plot(tspan,p_traj.x_0 + p_traj.r * cos(p_traj.omega*tspan) ,'r--');
    plot(tspan,rE1(2,:),'b','LineWidth',2)
    % plot(tspan,p_traj.y_0 + p_traj.r * sin(p_traj.omega*tspan) ,'b--');
    plot(tspan,rE2(1,:),'m','LineWidth',2)
    plot(tspan,rE2(2,:),'c','LineWidth',2)
    hold off %added 

    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x1','x2','y1','y2'});
    % xlabel('Time (s)'); ylabel('Position (m)'); legend({'x1','x_d','y1','y_d'});
    figure(3); clf;
    plot(tspan,vE1(1,:),'r','LineWidth',2)
    hold on
    plot(tspan,vE1(2,:),'b','LineWidth',2)
    hold on
    plot(tspan,vE2(1,:),'m','LineWidth',2)
    hold on
    plot(tspan,vE2(2,:),'c','LineWidth',2)
    
    xlabel('Time (s)'); ylabel('Velocity (m)'); legend({'vel_x1','vel_y1','vel_x2','vel_y2'});
    
    figure(4)
    plot(tspan,z_out(1:4,:)*180/pi) %convert gen coord: rad to degrees
    legend('q1','q2','q3','q4');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    
    figure(5)
    plot(tspan,z_out(5:8,:)*180/pi) %convert gen vel's: rad/s to degrees/s
    legend('q1dot','q2dot','q3dot','q4dot');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/sec)');
    
    %% Animate Solution
    figure(6); clf;
    hold on
  
    % Target traj
    TH = 0:.1:2*pi;
    plot( p_traj.x_0 + p_traj.r * cos(TH), ...
          p_traj.y_0 + p_traj.r * sin(TH),'k--'); 
    
    % Ground Q2.3
    plot([-.2 .2],[ground_height ground_height],'k'); 
    
    animateSol(tspan, z_out,p);
end

function tau = control_law(t, z, p, p_traj)
    % Controller gains, Update as necessary for Problem 1
    K_x = 150.; % Spring stiffness X
    K_y = 150.; % Spring stiffness Y
    D_x = 10.;  % Damping X
    D_y = 10.;  % Damping Y

    K = [K_x, 0 ; 0, K_y]; %K_xy = 0 
    D = [D_x, 0 ; 0, D_y]; %D_xy = 0 

    % Desired position of foot is a circle
    omega_swing = p_traj.omega; %rad/sec
    rEd1 = [p_traj.x_0 p_traj.y_0 0]' + ...
            p_traj.r*[cos(omega_swing*t) sin(omega_swing*t) 0]';
    % Compute desired velocity of foot
    vEd1 = p_traj.r*[-sin(omega_swing*t)*omega_swing    ...
                     cos(omega_swing*t)*omega_swing   0]';
    % Desired acceleration
    aEd1 = p_traj.r*[-cos(omega_swing*t)*omega_swing^2 ...
                    -sin(omega_swing*t)*omega_swing^2 0]';
    % Leg 2 
    PS = pi/2; % Phase shift for legs (rad)! legs move at 90deg out-of-phase 
    rEd2 = [p_traj.x_0 p_traj.y_0 0]' + ...
            p_traj.r*[cos(omega_swing*(t+PS)) sin(omega_swing*(t+PS)) 0]';
    vEd2 = p_traj.r*[-sin(omega_swing*(t+PS))*omega_swing    ...
                     cos(omega_swing*(t+PS))*omega_swing   0]';
    aEd2 = p_traj.r*[-cos(omega_swing*(t+PS))*omega_swing^2 ...
                    -sin(omega_swing*(t+PS))*omega_swing^2 0]';
    
    % Actual position and velocity 
    rE1 = position_foot1(z,p);
    vE1 = velocity_foot1(z,p);
    rE2 = position_foot2(z,p);
    vE2 = velocity_foot2(z,p);

    % To get f, WE NEED LAMBDA, MU, & RHO
    lambda = inv (jacobian_foot1(z,p) * inv(A_leg(z,p)) * (jacobian_foot1(z,p).'));
    mu1 = (lambda * jacobian_foot1(z,p) * inv(A_leg(z,p)) * Corr_leg(z,p)) - (lambda * jacobian_dot_foot1(z,p) * qdot_twolegs(z,p));
    rho1 = lambda * jacobian_foot1(z,p) * inv(A_leg(z,p)) * Grav_leg(z,p);
    % Compute virtual foce 
    intermediate1 = aEd1 + K * (rEd1 - rE1(1:2, 1:end)) + D * (vEd1 - vE1(1:2, 1:end));
    % intermediate2 = K * (rEd1 - rE1(1:2, 1:end)) + D * (vEd1 - vE1(1:2, 1:end)); %without aceeleration component
    f1 = lambda * intermediate1 + mu1 + rho1; % should result in a 2x1 matrix

    lambda2 = inv (jacobian_foot2(z,p) * inv(A_leg(z,p)) * (jacobian_foot2(z,p).'));
    mu2 = (lambda2 * jacobian_foot2(z,p) * inv(A_leg(z,p)) * Corr_leg(z,p)) - (lambda2 * jacobian_dot_foot2(z,p) * qdot_twolegs(z,p));
    rho2 = lambda2 * jacobian_foot2(z,p) * inv(A_leg(z,p)) * Grav_leg(z,p);
    intermediate2 = aEd2 + K * (rEd2 - rE2(1:2, 1:end)) + D * (vEd2 - vE2(1:2, 1:end));
    f2 = lambda2 * intermediate2 + mu2 + rho2;

    % % Compute virtual foce for Question 1.4 and 1.5
    % f  = [K_x * (rEd1(1) - rE(1) ) + D_x * ( - vE(1) ) ;
    %       K_y * (rEd1(2) - rE(2) ) + D_y * ( - vE(2) ) ];
    % 
    %% Task-space compensation and feed forward for Question 1.8

    % Map to joint torques  
    J1  = jacobian_foot1(z,p);
    tau1 = J1' * f1;
    J2  = jacobian_foot2(z,p);
    tau2 = J2' * f2;
    tau = [tau1, tau2];
end


function dz = dynamics(t,z,p,p_traj)
    % Get mass matrix
    A = A_twolegs(z,p);
    
    % Compute Controls
    tau = control_law(t,z,p,p_traj);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_twolegs(z,tau,p);
    
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1:4) = z(5:8);
    dz(5:8) = qdd;
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)
    % Actual position and velocity 
    rE1 = position_foot1(z,p);
    vE1 = velocity_foot1(z,p);
    rE2 = position_foot2(z,p);
    vE2 = velocity_foot2(z,p);
    J1  = jacobian_foot1(z,p); % 2x2 matrix
    J2  = jacobian_foot2(z,p); % 2x2 matrix
    zero = zeros(2,2);
    J = [J1, zero; zero, J2];%4x4 matrix
    osim1 = inv (J * inv(A_twolegs(z,p)) * J); 
    mass_eff_x1 = 1 / ([1,0,0,0] * inv(osim1) * [1;0;0;0]); %m_eff_r = 1/(rhat*inv(osim)*rhat.')
    mass_eff_y1 = 1 / ([0,1,0,0] * inv(osim1) * [0;1;0;0]);

    C_y1 = rE1(2,1:end) - yC;
    %Cdot_y = vE(2,1:end);
    Cdot_y1 = vE1(2);

    qdot = qdot_twolegs(z,p);

    if C_y1 < 0 & Cdot_y1 < 0 
        disp('entered')
        F_impulse_y = mass_eff_y1 * (-rest_coeff * Cdot_y1 - J1(2,:) * qdot_twolegs(z,p));
        %first qdot update --> updated later 
        %qdot = qdot + (inv(A_leg(z,p)) * J(2,:).' + F_impulse_y);
        F_impulse_x = mass_eff_x1 * (-J(1,:) * qdot_twolegs(z,p));
        if abs(F_impulse_x) > abs(fric_coeff * F_impulse_y)
            F_impulse_x = sign(F_impulse_x) * fric_coeff * F_impulse_y;
        end
        %update both for x and y impulses
        qdot = qdot + (inv(A_leg(z,p)) * J(2,:).' * F_impulse_y) + (inv(A_leg(z,p)) * J(1,:).' * F_impulse_x);
        % qdot = qdot_leg(z,p);
    end
    qdot = qdot;

end

function qdot = joint_limit_constraint(z,p)

end

function animateSol(tspan, x,p)
    % Prepare plot handles
    hold on
    h_OB = plot([0],[0],'LineWidth',2);
    h_AC = plot([0],[0],'LineWidth',2);
    h_BD = plot([0],[0],'LineWidth',2);
    h_CE = plot([0],[0],'LineWidth',2);
   
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-.2 .2 -.3 .1]);

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end
        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_leg(z,p);

        rA = keypoints(:,1); % Vector to base of cart
        rB = keypoints(:,2);
        rC = keypoints(:,3); % Vector to tip of pendulum
        rD = keypoints(:,4);
        rE = keypoints(:,5);

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        set(h_OB,'XData',[0 rB(1)]);
        set(h_OB,'YData',[0 rB(2)]);
        
        set(h_AC,'XData',[rA(1) rC(1)]);
        set(h_AC,'YData',[rA(2) rC(2)]);
        
        set(h_BD,'XData',[rB(1) rD(1)]);
        set(h_BD,'YData',[rB(2) rD(2)]);
        
        set(h_CE,'XData',[rC(1) rE(1)]);
        set(h_CE,'YData',[rC(2) rE(2)]);

        pause(.01)
    end
end