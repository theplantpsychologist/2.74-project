syms L
g = 9.81;
Patm = 1e5;
pw = 1000;
pm = 13560;
hw2 = 18e-2;
hm1 = 32e-2;
hm2 = 15e-2;
hw = L*cos(35/180*pi);
P1 = Patm + pw*(hw-hw2)*g;
P2 = P1 + pm*(hm1-hm2)*g;
Pa = subs(P2,L,120e-2);
% eq = P2 == 135e3;
% solve(eq, L)
Pa

%%
clc
clear
syms h
g = 9.81;
Pair = 110e3;
x1 = 8e-2;
Pair2 = 150e3;
x2 = Pair*x1/Pair2
xwater = 9e-2;
xtotal = x1+xwater+0.12;
xm = xtotal-xwater-x2
pw = 998;
pm = 13550;
T = 293;
P2 = Pair2+(9e-2)*pw*g;
P4 = P2 + xm*pm*g;
Patm = 101e3;
eq = Patm == P4-pm*g*h;
solve(eq,h)
% (Pair+9e-2*pw*g+0.12*pm*g-Patm)/(pm*g);

%%
clc
clear
zmax = 15;
pw = 1030;
w = 5;
g = 9.81;
sin_theta = 0.6;
zeta_max = 10
F = (pw*g*(zmax-9)^2/(2*sin_theta)+9*((zmax-9)/sin_theta))*w
F2 = (pw*g*sin_theta*zeta_max^2/2+9*zeta_max)*w