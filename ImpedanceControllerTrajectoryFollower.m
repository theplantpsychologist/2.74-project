% This is the main MATLAB script for Lab 5.
%
% You will need to modify the Mbed code and this script, but should not need to make any other changes.
%
%% SET YOUR INPUTS HERE

% Bezier curve control points
const_point = [0; -0.15]; %[x;y] or [q1,q2] constant coordinate (x,q1,q2 coordinates should be opposite sign due to direction motors are mounted)
pts_foot = repmat(const_point,1,8);

% Circle parameters
radius = 0.07;             % Radius of the circle (e.g., 5 cm)
x_center = 0;        % X-coordinate of the center
y_center = -0.15;        % Y-coordinate of the center
num_points = 8;            % Number of points to approximate the circle

% Generate angles for circular trajectory
angles = linspace(0, 2*pi, num_points + 1);  % 0 to 2π, with num_points divisions
angles(end) = [];  % Remove the last point to avoid overlap at the start

% Calculate circular trajectory points
x_pts = x_center + radius * cos(angles);  % X-coordinates
y_pts = y_center + radius * sin(angles);  % Y-coordinates

% Combine into pts_foot matrix
pts_foot_1 = [x_pts; y_pts];
pts_foot_1 = [x_pts x_pts; y_pts y_pts];
        
% Initial leg angles for encoder resets (negative of q1,q2 in lab handout due to direction motors are mounted)
angle1_init = 0;
angle2_init = -pi/2; 

% Total experiment time is buffer,trajectory,buffer
traj_time         = 5;
pre_buffer_time   = 2; % this should be 0 for constant points, 2 for Bezier trajectories
post_buffer_time  = 1;

% Gains for impedance controller
% If a gain is not being used in your Mbed code, set it to zero
% For joint space control, use K_xx for K1, K_yy for K2, D_xx for D1, D_yy for D2
gains.K_xx = 300;
gains.K_yy = 300;
gains.K_xy = 0;

gains.D_xx = 5;
gains.D_yy = 5;
gains.D_xy = 0;

% Maximum duty cycle commanded by controller (should always be <=1.0)
duty_max   = 0.4;

%% Run Experiment
[output_data] = RunTrajectoryExperiment(angle1_init, angle2_init, pts_foot_1,...
                                        traj_time, pre_buffer_time, post_buffer_time,...
                                        gains, duty_max);

%% Extract data
t = output_data(:,1);
x_1 = -output_data(:,12); % actual foot position in X (negative due to direction motors are mounted)
y_1 = output_data(:,13); % actual foot position in Y
   
xdes_1 = -output_data(:,16); % desired foot position in X (negative due to direction motors are mounted)
ydes_1 = output_data(:,17); % desired foot position in Y

%% Plot foot vs desired
figure(3); clf;
subplot(211); hold on
plot(t,xdes_1,'r-'); plot(t,x_1);
xlabel('Time (s)'); ylabel('X (m)'); legend({'Desired','Actual'});

subplot(212); hold on
plot(t,ydes_1,'r-'); plot(t,y_1);
xlabel('Time (s)'); ylabel('Y (m)'); legend({'Desired','Actual'});

figure(4); clf; hold on
plot(xdes_1,ydes_1,'r-'); plot(x_1,y_1,'k');
xlabel('X (m)'); ylabel('Y (m)'); legend({'Desired','Actual'});
