% h = 200;
% l = 25e-3;
% d_out = 5e-3;
% t = 1e-3;
% k = 0.2;
% d_in = d_out - 2*t;
% V = pi*d_out^2/4*l -pi*d_in^2/4*l
% A = d_out*pi*l + d_in*pi*l
% Lc = V/A
% Bi = h*Lc/k

cyl_V = @(r,l) r^2*pi*l;
cyl_A = @(r,l) r*2*pi*l;
sph_V = @(r)  4/3*r^3*pi;
sph_A = @(r)  4*pi*r^2;

% k = 318;
% h = 20;
% d = 33e-3;
% t = 1e-3;
% A = (d/2)^2*pi*2
% V = A/2*t

% k = 3;
% h = 5;
% r = 58e-2/2;
% A = sph_A(r);
% V = sph_V(r);

% k = 16;
% h = 200;
% r = 1e-2/2;
% l = 15e-2;
% A = cyl_A(r,l);
% V = cyl_V(r,l);

% k = 0.3;
% ha = 10;
% hw = 100;
% r = 10;
% l = 100;
% A = cyl_A(r,l)
% V = cyl_V(r,l)

% h = 50
% k = 235
% V = 100e-6
% A = 1000e-4
% Lc = V/A
% Bi = Lc*h/k
% 
% 
% T0 = 15 + 273;
% Ts = -10 + 273;
% T_target = 273;
% k = 0.52;
% rho = 2050;
% c = 1840; 
% x = 1;
% alpha = k/(rho*c);
% t = x^2/(4*alpha*erfcinv((T_target-Ts)/(T0-Ts))^2)

%%
clear
cyl_V = @(r,l) r^2*pi*l;
cyl_A = @(r,l) r*2*pi*l;
sph_V = @(r)  4/3*r^3*pi;
sph_A = @(r)  4*pi*r^2;

% k = 400;
% c = 384;
% rho = 8940;
% R = 77e-6;
% d = 0.01;
% l = 1;
% I = 1000;
% Lc = d/4;
% h = 300;
% Bi = h*Lc/k;
% V = cyl_V(d/2,l);
% A = cyl_A(d/2,l);
% tau = rho*c*V/h/A;
% Tinf = 20+273;
% Qgen = I^2*R;
% Teq = Qgen/h/A+Tinf;
% 
% m = rho*V
% S = m*c*log(Tinf/Teq)
Ti = 601;
Tinf = 300;
h = 240;
rho = 11340;
c = 129;
k = 34;
d = 3e-3;
r = d/2;
V = sph_V(r);
A = sph_A(r);
Lc = V/A;
Bi = h*Lc/k;
tau = rho*V*c/h/A;
theta = Ti-Tinf;
T_target = 450;
theta_0 = T_target-Tinf;
t = log(theta/theta_0)*tau;
g = 9.81;
x = 1/2*g*t^2;
m = rho*V;
Sair = m*c*log(Tinf/Ti)
Q = m*c*(Ti-Tinf);
Ssph = Q/Tinf
Sgen = Sair-Ssph

%%
clear
k = 0.17;
alpha = 1.28e-7;
Ti = 303;
Ts = 623;
Tign = 583;
t = 14400;
h = 65;
Lc = 0.2;
Bi = h*Lc/k;
syms T
theta = T - Ts;
theta_0 = Ti-Ts;
n = 1;
lambda = (n+1/2)*pi;
eq = theta/theta_0 == 4*sin(lambda)^2/(lambda*(2*lambda+sin(2*lambda)))*exp(-lambda^2*alpha*t/(Lc^2));
vpa(solve(eq, T),4)
%%
clear
T = 50+273;
k_m = 13;
rho_m = 8000;
c_m = 470;

Ti = 37;
k_t = 0.6;
rho_t = 1000;
c_t = 4800;
Tdeath = 48+273;
tdeath = 10; 

alpha_t = k_t/rho_t/c_t

%%
clear 
% Given data for tissue and machinery properties
k_tissue = 0.6;             % W/m-K
rho_tissue = 1000;          % kg/m^3
c_tissue = 4800;            % J/kg-K
alpha_tissue = k_tissue / (rho_tissue * c_tissue);  % m^2/s

k_machine = 13;             % W/m-K
rho_machine = 8000;         % kg/m^3
c_machine = 470;            % J/kg-K
T_i = 37;                   % Initial temperature of tissue in °C
T_target = 48;              % Target temperature for injury in °C

% Machinery temperature and time ranges with finer steps
T_machine = 50:1:100;       % °C, finer machinery temperature range
time_range = 10:0.5:30;     % seconds, finer time range for contact duration

% Create a mesh grid for T_machine and time_range
[T_m_grid, time_grid] = meshgrid(T_machine, time_range);

% Initialize depth of injury matrix
depth_injury = zeros(size(T_m_grid));

% Loop over each (T_machine, time) combination in the mesh grid
for i = 1:size(T_m_grid, 1)
    for j = 1:size(T_m_grid, 2)
        % Current machinery temperature and contact time
        T_m = T_m_grid(i, j);
        t = time_grid(i, j);
        
        % Surface temperature calculation
        T_s = (T_i * sqrt(k_tissue * rho_tissue * c_tissue) + T_m * sqrt(k_machine * rho_machine * c_machine)) ...
              / (sqrt(k_tissue * rho_tissue * c_tissue) + sqrt(k_machine * rho_machine * c_machine));
        
        % Solving for depth x where T(x, t) = T_target
        % Using inverse of erfc to find x
        theta = (T_target - T_i) / (T_s - T_i);
        if theta > 0 && theta < 1  % Only calculate for valid theta values
            z = erfcinv(theta);  % Solve for erfc argument
            x = 2 * z * sqrt(alpha_tissue * t);  % Depth x for the target temperature
            depth_injury(i, j) = x;  % Store depth of injury for this (T_m, t)
        else
            depth_injury(i, j) = NaN;  % Out of valid range for erfc, set as NaN
        end
    end
end

% Plot the results
figure;
mesh(T_m_grid, time_grid, depth_injury);
xlabel('Machinery Temperature (°C)');
ylabel('Contact Time (s)');
zlabel('Depth of Injury (m)');
title('Depth of Irreversible Tissue Damage (Finer Mesh)');
colorbar;

