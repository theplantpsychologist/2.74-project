function E = energy_twolegs(in1,in2)
%ENERGY_TWOLEGS
%    E = ENERGY_TWOLEGS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    04-Nov-2024 14:09:58

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I4 = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
dth1 = in1(5,:);
dth2 = in1(6,:);
dth3 = in1(7,:);
dth4 = in1(8,:);
g = in2(19,:);
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
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = cos(th1);
t3 = cos(th2);
t4 = cos(th3);
t5 = cos(th4);
t6 = th1+th2;
t7 = th3+th4;
t8 = N.^2;
t9 = dth1.^2;
t10 = dth2.^2;
t11 = dth3.^2;
t12 = dth4.^2;
t13 = l_AC.^2;
t14 = l_A_m3.^2;
t15 = l_B_m2.^2;
t16 = l_C_m4.^2;
t17 = l_OA.^2;
t18 = l_OB.^2;
t19 = l_O_m1.^2;
t20 = l_OA.*t2;
t21 = l_OA.*t4;
t22 = cos(t6);
t23 = cos(t7);
et1 = (I1.*t9)./2.0+(I2.*t9)./2.0+(I1.*t11)./2.0+(I2.*t10)./2.0+(I3.*t9)./2.0+(I2.*t11)./2.0+(I3.*t10)./2.0+(I4.*t9)./2.0+(I2.*t12)./2.0+(I3.*t11)./2.0+(I3.*t12)./2.0+(I4.*t11)./2.0+(Ir.*t9)./2.0+(Ir.*t11)./2.0-g.*m2.*(l_B_m2.*t22+l_OB.*t2)-g.*m2.*(l_B_m2.*t23+l_OB.*t4)-g.*m3.*(t20+l_A_m3.*t22)-g.*m3.*(t21+l_A_m3.*t23)+I2.*dth1.*dth2+I3.*dth1.*dth2+I2.*dth3.*dth4+I3.*dth3.*dth4-g.*m4.*(t20+l_AC.*t22+l_C_m4.*t2)-g.*m4.*(t21+l_AC.*t23+l_C_m4.*t4)+(Ir.*t8.*t9)./2.0+(Ir.*t8.*t10)./2.0+(Ir.*t8.*t11)./2.0+(Ir.*t8.*t12)./2.0;
et2 = (m2.*t9.*t15)./2.0+(m3.*t9.*t14)./2.0+(m4.*t9.*t13)./2.0+(m2.*t10.*t15)./2.0+(m3.*t10.*t14)./2.0+(m4.*t10.*t13)./2.0+(m2.*t11.*t15)./2.0+(m3.*t11.*t14)./2.0+(m4.*t11.*t13)./2.0+(m1.*t9.*t19)./2.0+(m2.*t9.*t18)./2.0+(m2.*t12.*t15)./2.0+(m3.*t9.*t17)./2.0+(m3.*t12.*t14)./2.0+(m4.*t9.*t16)./2.0+(m4.*t12.*t13)./2.0+(m4.*t9.*t17)./2.0+(m1.*t11.*t19)./2.0+(m2.*t11.*t18)./2.0+(m3.*t11.*t17)./2.0+(m4.*t11.*t16)./2.0+(m4.*t11.*t17)./2.0+dth1.*dth2.*m2.*t15+dth1.*dth2.*m3.*t14+dth1.*dth2.*m4.*t13+dth3.*dth4.*m2.*t15+dth3.*dth4.*m3.*t14+dth3.*dth4.*m4.*t13-g.*l_O_m1.*m1.*t2-g.*l_O_m1.*m1.*t4;
et3 = l_C_m4.*l_OA.*m4.*t9+l_C_m4.*l_OA.*m4.*t11+Ir.*N.*dth1.*dth2+Ir.*N.*dth3.*dth4+l_AC.*l_C_m4.*m4.*t3.*t9+l_AC.*l_C_m4.*m4.*t5.*t11+l_AC.*l_OA.*m4.*t3.*t9+l_AC.*l_OA.*m4.*t5.*t11+l_A_m3.*l_OA.*m3.*t3.*t9+l_A_m3.*l_OA.*m3.*t5.*t11+l_B_m2.*l_OB.*m2.*t3.*t9+l_B_m2.*l_OB.*m2.*t5.*t11+dth1.*dth2.*l_AC.*l_C_m4.*m4.*t3+dth3.*dth4.*l_AC.*l_C_m4.*m4.*t5+dth1.*dth2.*l_AC.*l_OA.*m4.*t3+dth3.*dth4.*l_AC.*l_OA.*m4.*t5+dth1.*dth2.*l_A_m3.*l_OA.*m3.*t3+dth3.*dth4.*l_A_m3.*l_OA.*m3.*t5+dth1.*dth2.*l_B_m2.*l_OB.*m2.*t3+dth3.*dth4.*l_B_m2.*l_OB.*m2.*t5;
E = et1+et2+et3;
end
