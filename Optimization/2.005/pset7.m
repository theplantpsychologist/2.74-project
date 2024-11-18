
Ta = 45;
Tb = 45;
Tc = 0;
Tinf = 25;
ha = 10;
hb = 100;
hc = 1000;
ka = 100;
kb = 50;
kc = 1;
L = 0.5;

qa = (Ta-Tinf)/(1/ha+L/ka)
qb = (Tb-Tinf)/(1/hb+L/kb)
qc = (Tc-Tinf)/(1/hc+L/kc)
Ts_a = qa/ha+Tinf
Ts_b = qb/hb+Tinf
Ts_c = qc/hc+Tinf
%%
% La = 0.3;
% Lb = 0.1;
% Lc = 0.2;
% kb = 100;
% kc = 10;
% ka = 25;
% 
% R_ = La/ka + Lb/kb + Lc/kc;
% a = ka/R_
% b = kb/R_

% Ac = 1e-4;
% P = 0.04;
% h = 10;
% k = 10;
% L = 0.05;
% m = sqrt(h*P/(k*Ac))

% 
% Ac = 2.5e-5;
% P = 0.02;
% h = 100;
% k = 10;
% L = 0.05;
% m = sqrt(h*P/(k*Ac))
% 
% Ac = 1e-4;
% P = 0.04;
% h = 1000;
% k = 10;
% L = 0.01;
% m = sqrt(h*P/(k*Ac))
% 
% m*L
% h/(m*k)

clc
clear
xA= 0.15;
xB = 0.075;
kA = 70;
kb = (xB^2)/(xA^2)*kA
%%
clc
clear
Tc = 9.2;
kcore = 0.847;
kn = 67.3;
ks = 18.2;
kcu = 1414;
Tl = 1.3;
Th = 17;
l = 25e-2;
Acore = 7.854e-9;
Ashell = 1.62e-8;
Rmg = l/(kcore*Acore);

syms dQ Rtotal Ln
Rn = (l-Ln)/(ks*Ashell)+Ln/(kn*Ashell);
Rnorm = 1/(Acore*kcore/Ln+ kn*Ashell/Ln);
eq1 = dQ == (Th - Tl)/Rtotal;
eq2 = Rtotal == 1/(1/Rmg + 1/Rn);
% eq3 = dQ == (Th - Tc)/Ln*(Acore*kcore+kn*Ashell);
eq3 = dQ == (Th-Tc)/Rnorm;
[dQ_ans, Rtotal_ans, Ln_ans] = solve([eq1, eq2, eq3], [dQ, Rtotal,Ln]);
dQ_ans = vpa(dQ_ans, 4)
Rtotal_ans = vpa(Rtotal_ans, 4)
Ln_ans = vpa(Ln_ans, 4)

%%
clc
clear
%parameters

tal_in = 1e-3;
tal_out = 2e-3;
tperl = 1.5e-2;

rn = 15e-2;
r1 = rn+tal_in;
r2 = r1+tperl;
r3 = r2+tal_out;

Tair = 300;
Tn = 77;

kal = 230;
kperl = 0.04;
hair = 5;
hr = 1;
hn = 1000;
L = 198e3;

%find areas
vol = @(r) 4*pi*r^2;
An = vol(rn);
Aal_out = vol(r3);

%find thermal resistances
Rconv_n = 1/hn*An;

Rcond = @(r_out,r_in,k) (1/r_in+1/r_out)/(4*pi*k);
Rcond_al_in = Rcond(r1, rn, kal);
Rcond_perl = Rcond(r2, r1, kperl);
Rcond_al_out = Rcond(r3, r2, kal);

Rconv_air = 1/(hair*Aal_out);
Rrad_air = 1/(hr*Aal_out);
Rair = 1/(1/Rconv_air+1/Rrad_air);

Rtotal = Rconv_n+Rcond_al_in+Rcond_perl+Rcond_al_out+Rair

%find heat transfer rate
dQ = (Tn-Tair)/Rtotal

%find evaporation rate
dm = dQ/L
%%
dQ = 20*2*600e-4
dQ*20e6*3.154e7*550/(3.6e6)
