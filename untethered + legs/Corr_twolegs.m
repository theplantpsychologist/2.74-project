function Corr_Joint_Sp = Corr_twolegs(in1,in2)
%Corr_twolegs
%    Corr_Joint_Sp = Corr_twolegs(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    20-Nov-2024 15:05:55

dth1 = in1(9,:);
dth2 = in1(10,:);
dth3 = in1(11,:);
dth4 = in1(12,:);
l_AC = in2(17,:);
l_A_m3 = in2(13,:);
l_B_m2 = in2(12,:);
l_C_m4 = in2(14,:);
l_OA = in2(15,:);
l_OB = in2(16,:);
l_O_m1 = in2(11,:);
m1 = in2(1,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
th4 = in1(6,:);
t2 = cos(th1);
t3 = cos(th3);
t4 = sin(th1);
t5 = sin(th2);
t6 = sin(th3);
t7 = sin(th4);
t8 = th1+th2;
t9 = th3+th4;
t10 = dth1.^2;
t11 = dth2.^2;
t12 = dth3.^2;
t13 = dth4.^2;
t14 = l_AC.*l_C_m4.*m4;
t15 = l_AC.*l_OA.*m4;
t16 = l_A_m3.*l_OA.*m3;
t17 = l_B_m2.*l_OB.*m2;
t18 = cos(t8);
t19 = cos(t9);
t20 = sin(t8);
t21 = sin(t9);
t22 = t14+t15+t16+t17;
et1 = -l_AC.*m4.*t10.*t20-l_AC.*m4.*t11.*t20-l_AC.*m4.*t12.*t21-l_AC.*m4.*t13.*t21-l_A_m3.*m3.*t10.*t20-l_A_m3.*m3.*t11.*t20-l_A_m3.*m3.*t12.*t21-l_A_m3.*m3.*t13.*t21-l_B_m2.*m2.*t10.*t20-l_B_m2.*m2.*t11.*t20-l_B_m2.*m2.*t12.*t21-l_B_m2.*m2.*t13.*t21-l_C_m4.*m4.*t4.*t10-l_C_m4.*m4.*t6.*t12-l_OA.*m3.*t4.*t10-l_OB.*m2.*t4.*t10-l_OA.*m4.*t4.*t10-l_OA.*m3.*t6.*t12-l_OB.*m2.*t6.*t12-l_OA.*m4.*t6.*t12-l_O_m1.*m1.*t4.*t10-l_O_m1.*m1.*t6.*t12-dth1.*dth2.*l_AC.*m4.*t20.*2.0-dth3.*dth4.*l_AC.*m4.*t21.*2.0-dth1.*dth2.*l_A_m3.*m3.*t20.*2.0-dth3.*dth4.*l_A_m3.*m3.*t21.*2.0-dth1.*dth2.*l_B_m2.*m2.*t20.*2.0;
et2 = dth3.*dth4.*l_B_m2.*m2.*t21.*-2.0;
Corr_Joint_Sp = [et1+et2;l_AC.*m4.*t10.*t18+l_AC.*m4.*t11.*t18+l_AC.*m4.*t12.*t19+l_AC.*m4.*t13.*t19+l_A_m3.*m3.*t10.*t18+l_A_m3.*m3.*t11.*t18+l_A_m3.*m3.*t12.*t19+l_A_m3.*m3.*t13.*t19+l_B_m2.*m2.*t10.*t18+l_B_m2.*m2.*t11.*t18+l_B_m2.*m2.*t12.*t19+l_B_m2.*m2.*t13.*t19+l_C_m4.*m4.*t2.*t10+l_C_m4.*m4.*t3.*t12+l_OA.*m3.*t2.*t10+l_OB.*m2.*t2.*t10+l_OA.*m4.*t2.*t10+l_OA.*m3.*t3.*t12+l_OB.*m2.*t3.*t12+l_OA.*m4.*t3.*t12+l_O_m1.*m1.*t2.*t10+l_O_m1.*m1.*t3.*t12+dth1.*dth2.*l_AC.*m4.*t18.*2.0+dth3.*dth4.*l_AC.*m4.*t19.*2.0+dth1.*dth2.*l_A_m3.*m3.*t18.*2.0+dth3.*dth4.*l_A_m3.*m3.*t19.*2.0+dth1.*dth2.*l_B_m2.*m2.*t18.*2.0+dth3.*dth4.*l_B_m2.*m2.*t19.*2.0;-dth2.*t5.*t22.*(dth1.*2.0+dth2);t5.*t10.*t22;-dth4.*t7.*t22.*(dth3.*2.0+dth4);t7.*t12.*t22];
end
