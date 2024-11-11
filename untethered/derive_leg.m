clear
name = 'leg';

% Define variables for time, generalized coordinates + derivatives, controls, and parameters 
syms t x y dx dy ddx ddy th1 th2 dth1 dth2 ddth1 ddth2 real
syms m1 m2 m3 m4 I1 I2 I3 I4 l_O_m1 l_B_m2 l_A_m3 l_C_m4  g real
syms l_OA l_OB l_AC l_DE real 
syms tau1 tau2 Fx Fy real
syms Ir N real

% Group them
q   = [x; y; th1  ; th2 ];      % generalized coordinates
dq  = [dx; dy; dth1 ; dth2];    % first time derivatives
ddq = [ddx; ddy; ddth1; ddth2];  % second time derivatives
u   = [tau1 ; tau2; Fx ; Fy];     % controls


p   = [m1 m2 m3 m4 I1 I2 I3 I4 Ir N l_O_m1 l_B_m2 l_A_m3 l_C_m4 l_OA l_OB l_AC l_DE g]';        % parameters

% Generate Vectors and Derivativess
ihat = [0; -1; 0];
jhat = [1; 0; 0];

xhat = [1; 0; 0];
yhat = [0; 1; 0];

khat = cross(ihat,jhat);
e1hat =  cos(th1)*ihat + sin(th1)*jhat;
e2hat =  cos(th1+th2)*ihat + sin(th1+th2)*jhat;

ddt = @(r) jacobian(r,[q;dq])*[dq;ddq]; % a handy anonymous function for taking time derivatives

rO = x*xhat+y*yhat;
rA = rO + l_OA * e1hat;
rB = rO + l_OB * e1hat;
rC = rA  + l_AC * e2hat;
rD = rB  + l_AC * e2hat;
rE = rD  + l_DE * e1hat;

r_m1 = rO + l_O_m1 * e1hat;
r_m2 = rO + rB + l_B_m2 * e2hat;
r_m3 = rO + rA + l_A_m3 * e2hat;
r_m4 = rO + rC + l_C_m4 * e1hat;

drA = ddt(rA);
drB = ddt(rB);
drC = ddt(rC);
drD = ddt(rD);
drE = ddt(rE);

dr_m1 = ddt(r_m1);
dr_m2 = ddt(r_m2);
dr_m3 = ddt(r_m3);
dr_m4 = ddt(r_m4);

% Calculate Kinetic Energy, Potential Energy, and Generalized Forces
F2Q = @(F,r) simplify(jacobian(r,q)'*(F));    % force contributions to generalized forces
M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));   % moment contributions to generalized forces


omega1 = dth1;
omega2 = dth1 + dth2;
omega3 = dth1 + dth2;
omega4 = dth1;

T1 = (1/2)*m1 * dot(dr_m1,dr_m1) + (1/2) * I1 * omega1^2;
T2 = (1/2)*m2 * dot(dr_m2,dr_m2) + (1/2) * I2 * omega2^2;
T3 = (1/2)*m3 * dot(dr_m3,dr_m3) + (1/2) * I3 * omega3^2;
T4 = (1/2)*m4 * dot(dr_m4,dr_m4) + (1/2) * I4 * omega4^2;
T1r = (1/2)*Ir*(N*dth1)^2;
T2r = (1/2)*Ir*(dth1 + N*dth2)^2;

Vg1 = m1*g*dot(r_m1, -ihat);
Vg2 = m2*g*dot(r_m2, -ihat);
Vg3 = m3*g*dot(r_m3, -ihat);
Vg4 = m4*g*dot(r_m4, -ihat);

T = simplify(T1 + T2 + T3 + T4 + T1r + T2r);
Vg = Vg1 + Vg2 + Vg3 + Vg4;

Q_tau1 = M2Q(tau1*khat,omega1*khat);
Q_tau2 = M2Q(tau2*khat,omega2*khat); 
Q_tau2R= M2Q(-tau2*khat,omega1*khat);
Q_tau = Q_tau1+Q_tau2 + Q_tau2R;

Q_force_x = F2Q(Fx*xhat, rO.*xhat);  % Forces acting at the base
Q_force_y = F2Q(Fy*yhat, rO.*yhat);
Q_force = Q_force_x + Q_force_y;

Q = Q_tau + Q_force;

% Assemble the array of cartesian coordinates of the key points
keypoints = [rO(1:2) rA(1:2) rB(1:2) rC(1:2) rD(1:2) rE(1:2)];

%% All the work is done!  Just turn the crank...
% Derive Energy Function and Equations of Motion
E = T+Vg;
L = T-Vg;
eom = ddt(jacobian(L,dq).') - jacobian(L,q).' - Q;



% Rearrange Equations of Motion
% A = simplify(jacobian(eom,ddq));
% b = A*ddq - eom;

% Modified code using equationsToMatrix
[A, b] = equationsToMatrix(eom, ddq);

% Equations of motion are
% eom = A *ddq + (coriolis term) + (gravitational term) - Q = 0
Mass_Joint_Sp = A;
Grav_Joint_Sp = simplify(jacobian(Vg, q)');
Corr_Joint_Sp = simplify( eom + Q - Grav_Joint_Sp - A*ddq);

% Compute foot jacobian
J = jacobian(rE,q);

% Compute ddt( J )
dJ= reshape( ddt(J(:)) , size(J) );

% Write Energy Function and Equations of Motion
z  = [q ; dq];

rE = rE(1:2);
drE= drE(1:2);
J  = J(1:2,1:4);
dJ = dJ(1:2,1:4);

matlabFunction(A,'file',['A_' name],'vars',{z p});
matlabFunction(b,'file',['b_' name],'vars',{z u p});
matlabFunction(E,'file',['energy_' name],'vars',{z p});
matlabFunction(rE,'file',['position_foot'],'vars',{z p});
matlabFunction(drE,'file',['velocity_foot'],'vars',{z p});
matlabFunction(J ,'file',['jacobian_foot'],'vars',{z p});
matlabFunction(dJ ,'file',['jacobian_dot_foot'],'vars',{z p});

matlabFunction(Grav_Joint_Sp ,'file', ['Grav_leg'] ,'vars',{z p});
matlabFunction(Corr_Joint_Sp ,'file', ['Corr_leg']     ,'vars',{z p});
matlabFunction(keypoints,'file',['keypoints_' name],'vars',{z p});

