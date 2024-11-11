function simulate_leg()
    %% Define fixed parameters
    m1 = .0393 + .2;      m2 = .0368; 
    m3 = .00783;          m4 = .0155;
    I1 = 25.1e-6;         I2 = 53.5e-6;
    I3 = 9.25e-6;         I4 = 22.176e-6;
    l_OA = .011;          l_OB = .042; 
    l_AC = .096;          l_DE = .091;
    l_O_m1 = 0.032;       l_B_m2 = 0.0344; 
    l_A_m3 = 0.0622;      l_C_m4 = 0.0610;
    N = 18.75;
    Ir = 0.0035 / N^2;
    g = 9.81;    

    q1_min = deg2rad(-40);
    q1_max = deg2rad(40);
    q2_min = deg2rad(20);
    q2_max = deg2rad(160);
    qlims = [q1_min, q1_max, q2_min, q2_max];

    restitution_coeff = 0;
    friction_coeff = 10;
    ground_height = -0.13;
    %% Parameter vector
    p = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 ...
         l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';
       
    %% Simulation Parameters
    dt = 0.001;
    tf = 2; % Adjusted for a shorter simulation time
    num_step = floor(tf / dt);
    tspan = linspace(0, tf, num_step); 
    % Initial state: [x; y; th1; th2; dx; dy; dth1; dth2]
    z0 = [0; l_OA+l_AC+l_DE; -pi / 4; pi / 2; 0; 0; 0; 0];
    z_out = zeros(8, num_step);
    z_out(:,1) = z0;

    for i = 1:num_step - 1
        dz = dynamics(tspan(i), z_out(:,i), p);
        % Update state with dynamics
        z_temp = z_out(:,i) + dz * dt;
    
        % Check for ground contact and update velocities if necessary
        z_temp(5:8) = discrete_impact_contact(z_temp, p, ...
            restitution_coeff, friction_coeff, ground_height);
    
        % Check for joint limit constraint and update velocities if necessary
        z_temp(5:8) = joint_limit_constraint(z_temp, p, qlims);
    
        % Update the state
        z_out(:,i+1) = z_temp;
    end

    %% Compute Energy
    E = energy_leg(z_out, p);
    figure(1); clf
    plot(tspan, E); xlabel('Time (s)'); ylabel('Energy (J)');

    %% Compute foot position over time
    rE = zeros(2, length(tspan));
    vE = zeros(2, length(tspan));
    for i = 1:length(tspan)
        rE(:, i) = position_foot(z_out(:, i), p);
        vE(:, i) = velocity_foot(z_out(:, i), p);
    end

    figure(2); clf;
    plot(tspan, rE(1,:), 'r', 'LineWidth', 2)
    hold on
    plot(tspan, rE(2,:), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Position (m)'); legend({'x', 'y'});

    figure(3); clf;
    plot(tspan, vE(1,:), 'r', 'LineWidth', 2)
    hold on
    plot(tspan, vE(2,:), 'b', 'LineWidth', 2)
    xlabel('Time (s)'); ylabel('Velocity (m/s)'); legend({'vel_x', 'vel_y'});

    figure(4)
    plot(tspan, z_out(3:4,:) * 180 / pi)
    legend('th1', 'th2');
    xlabel('Time (s)');
    ylabel('Angle (deg)');

    figure(5)
    plot(tspan, z_out(7:8,:) * 180 / pi)
    legend('dth1', 'dth2');
    xlabel('Time (s)');
    ylabel('Angular Velocity (deg/s)');

    %% Animate Solution
    figure(6); clf;
    hold on

    % Ground line
    plot([-.5 .5], [ground_height ground_height], 'k'); 

    animateSol(tspan, z_out, p);
end
function dz = dynamics(t, z, p)
    % Get mass matrix
    A = A_leg(z, p);

    % Control torques are zero (no actuation)
    tau = zeros(4, 1); % For [th1; th2] only

    % Get b = Q - V(q, dq) - G(q)
    b = b_leg(z, tau, p);

    % Solve for qdd.
    qdd = A \ b;
    dz = zeros(size(z));

    % Form dz
    dz(1:4) = z(5:8); % dq
    dz(5:8) = qdd;    % ddq
end

function qdot = discrete_impact_contact(z,p, rest_coeff, fric_coeff, yC)
    rE = position_foot(z,p);
    vE = velocity_foot(z,p);
    Cy = rE(2) - yC; %foot height relative to ground
    dCy = vE(2); %vertical speed
    qdot = z(5:8);% dx dy dth1 dth2

    if (Cy <= 0)&&(dCy<=0) %constraint violated
        J  = jacobian_foot(z,p);

        inv_mass = inv(A_leg(z,p));
        Jdot = jacobian_dot_foot(z,p);
        corr = Corr_leg(z,p); %corriolis effects
        grav = Grav_leg(z,p); %gravitational effects
        eff_mass = inv(J*inv_mass*J');
        ydir = [0;1];
        xdir = [1;0];
        osim_cy = ydir'*eff_mass*ydir; 
        osim_cx = xdir'*eff_mass*xdir;
        
        Jcy = J(2,:);
        Jcx = J(1,:);
        Fcy = osim_cy*(-rest_coeff*dCy - Jcy*qdot); % vertical ground reaction
        qdot = qdot + inv_mass*Jcy'*Fcy; %
        Fcx = osim_cx*(0-Jcx*qdot);
        
        if abs(Fcx) > fric_coeff*Fcy % if horizontal force is less than friction
            Fcx = fric_coeff*Fcy*Fcx/abs(Fcx); %truncate
        end
        qdot = qdot + inv_mass*Jcx'*Fcx;
        
    end

end
function qdot = joint_limit_constraint(z, p, limits)
    % Extract positions and velocities and limits
    q = z(1:4);      % Positions: [x; y; th1; th2]
    dq = z(5:8);     % Velocities: [dx; dy; dth1; dth2]
    qdot = dq;       % Initialize qdot with current velocities
    q1min = limits(1);
    q1max = limits(2);
    q2min = limits(3);
    q2max = limits(4);

    q1 = q(3);       % Joint angle th1
    q2 = q(4);       % Joint angle th2
    dth1 = dq(3);    % Joint velocity dth1
    dth2 = dq(4);    % Joint velocity dth2

    if ((q1 <= q1min) && (dth1 < 0)) || ((q1 >= q1max) && (dth1 > 0))% Joint limit violated and moving further into violation
        % Compute the Jacobian of the constraint
        % For joint limit on q1, the constraint is c(q) = q1 - q1_min >= 0
        % The gradient of c(q) with respect to q is:
        Jc = [0, 0, 1, 0]; % Only q1 has a non-zero partial derivative

        % Mass matrix
        A = A_leg(z, p);
        inv_mass = inv(A);

        % Compute the impulse needed to stop further violation
        % Constraint impulse: Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc')
        % For a perfectly inelastic collision (e = 0)
        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc');
        % -(1 + e) * (Jc * dq) is equivalent to -(gamma * dCy + Jc * dq)
        % because ddtCy = dCy/q*dq = Jc*dq
        
        % Adjust velocities
        qdot = qdot + inv_mass * (Jc' * Lambda);
    end

    if ((q2 <= q2min) && (dth2 < 0)) || ((q2 >= q2max) && (dth2 > 0))
        % c(q) = q2 - q2_max <= 0
        Jc = [0, 0, 0, 1];

        % Mass matrix
        A = A_leg(z, p);
        inv_mass = inv(A);

        e = 0; % Restitution coefficient for joint limit
        Lambda = - (1 + e) * (Jc * dq) / (Jc * inv_mass * Jc');

        % Adjust velocities
        qdot = qdot + inv_mass * (Jc' * Lambda);
    end
    
end


function animateSol(tspan, x, p)
    % Prepare plot handles
    hold on
    h_OB = plot([0], [0], 'LineWidth', 2);
    h_AC = plot([0], [0], 'LineWidth', 2);
    h_BD = plot([0], [0], 'LineWidth', 2);
    h_CE = plot([0], [0], 'LineWidth', 2); % CE segment

    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');

    axis equal
    axis([-.5 .5 -.3 .5]);

    % Step through and update animation
    for i = 1:length(tspan)
        % Skip frames to speed up the animation
        if mod(i, 10)
            continue;
        end
        t = tspan(i);
        z = x(:, i); 
        keypoints = keypoints_leg(z, p);

        rO = keypoints(:,1); % Vector to base of cart
        rA = keypoints(:,2);
        rB = keypoints(:,3); % Vector to tip of pendulum
        rC = keypoints(:,4);
        rD = keypoints(:,5);
        rE = keypoints(:,6);

        set(h_title, 'String', sprintf('t=%.2f', t)); % Update title


        % Update OB segment
        set(h_OB, 'XData', [rO(1), rB(1)]);
        set(h_OB, 'YData', [rO(2), rB(2)]);

        % Update AC segment
        set(h_AC, 'XData', [rA(1), rC(1)]);
        set(h_AC, 'YData', [rA(2), rC(2)]);

        % Update BD segment
        set(h_BD, 'XData', [rB(1), rD(1)]);
        set(h_BD, 'YData', [rB(2), rD(2)]);

        % Update CE segment
        set(h_CE, 'XData', [rC(1), rE(1)]);
        set(h_CE, 'YData', [rC(2), rE(2)]);

        pause(0.01)
    end
end

