function E = energy_twolegs(in1,in2)
%ENERGY_TWOLEGS
%    E = ENERGY_TWOLEGS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    23-Nov-2024 17:51:18

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I4 = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
dth1 = in1(9,:);
dth2 = in1(10,:);
dth3 = in1(11,:);
dth4 = in1(12,:);
dx = in1(7,:);
dy = in1(8,:);
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
th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
th4 = in1(6,:);
y = in1(2,:);
t2 = cos(th1);
t3 = cos(th3);
t4 = sin(th1);
t5 = sin(th3);
t6 = dth1+dth2;
t7 = dth3+dth4;
t8 = th1+th2;
t9 = th3+th4;
t10 = N.^2;
t11 = dth1.^2;
t12 = dth3.^2;
t25 = -y;
t13 = l_C_m4.*t2;
t14 = l_C_m4.*t3;
t15 = l_OA.*t2;
t16 = l_OB.*t2;
t17 = l_OA.*t3;
t18 = l_OB.*t3;
t19 = cos(t8);
t20 = cos(t9);
t21 = l_OA.*t4;
t22 = l_OA.*t5;
t23 = sin(t8);
t24 = sin(t9);
t26 = t6.^2;
t27 = t7.^2;
t28 = l_AC.*t19;
t29 = l_AC.*t20;
t30 = l_A_m3.*t19;
t31 = l_A_m3.*t20;
t32 = l_B_m2.*t19;
t33 = l_B_m2.*t20;
et1 = (I1.*t11)./2.0+(I1.*t12)./2.0+(I4.*t11)./2.0+(I4.*t12)./2.0+(I2.*t26)./2.0+(I2.*t27)./2.0+(I3.*t26)./2.0+(I3.*t27)./2.0+(m3.*((dx+dth1.*(t15+t30)+dth2.*t30).^2+(dy+dth1.*(t21+l_A_m3.*t23)+dth2.*l_A_m3.*t23).^2))./2.0+(m3.*((dx+dth3.*(t17+t31)+dth4.*t31).^2+(dy+dth3.*(t22+l_A_m3.*t24)+dth4.*l_A_m3.*t24).^2))./2.0+(m1.*((dx+dth1.*l_O_m1.*t2).^2+(dy+dth1.*l_O_m1.*t4).^2))./2.0+(m1.*((dx+dth3.*l_O_m1.*t3).^2+(dy+dth3.*l_O_m1.*t5).^2))./2.0+(Ir.*(dth1+N.*dth2).^2)./2.0+(Ir.*(dth3+N.*dth4).^2)./2.0;
et2 = (m4.*((dx+dth2.*t28+dth1.*(t13+t15+t28)).^2+(dy+dth1.*(t21+l_AC.*t23+l_C_m4.*t4)+dth2.*l_AC.*t23).^2))./2.0+(m4.*((dx+dth4.*t29+dth3.*(t14+t17+t29)).^2+(dy+dth3.*(t22+l_AC.*t24+l_C_m4.*t5)+dth4.*l_AC.*t24).^2))./2.0+(m2.*((dx+dth1.*(t16+t32)+dth2.*t32).^2+(dy+dth1.*(l_B_m2.*t23+l_OB.*t4)+dth2.*l_B_m2.*t23).^2))./2.0+(m2.*((dx+dth3.*(t18+t33)+dth4.*t33).^2+(dy+dth3.*(l_B_m2.*t24+l_OB.*t5)+dth4.*l_B_m2.*t24).^2))./2.0+g.*m1.*(y-l_O_m1.*t2)+g.*m1.*(y-l_O_m1.*t3)-g.*m4.*(t13+t15+t25+t28)-g.*m4.*(t14+t17+t25+t29)+(Ir.*t10.*t11)./2.0+(Ir.*t10.*t12)./2.0;
et3 = -g.*m3.*(t15+t25+t30)-g.*m2.*(t16+t25+t32)-g.*m3.*(t17+t25+t31)-g.*m2.*(t18+t25+t33);
E = et1+et2+et3;
end
