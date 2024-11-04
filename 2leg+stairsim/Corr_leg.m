function Corr_Joint_Sp = Corr_leg(in1,in2)
%Corr_leg
%    Corr_Joint_Sp = Corr_leg(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    04-Nov-2024 14:09:59

dth1 = in1(5,:);
dth2 = in1(6,:);
dth3 = in1(7,:);
dth4 = in1(8,:);
l_AC = in2(17,:);
l_A_m3 = in2(13,:);
l_B_m2 = in2(12,:);
l_C_m4 = in2(14,:);
l_OA = in2(15,:);
l_OB = in2(16,:);
m2 = in2(2,:);
m3 = in2(3,:);
m4 = in2(4,:);
th2 = in1(2,:);
th4 = in1(4,:);
t2 = sin(th2);
t3 = sin(th4);
t4 = l_AC.*l_C_m4.*m4;
t5 = l_AC.*l_OA.*m4;
t6 = l_A_m3.*l_OA.*m3;
t7 = l_B_m2.*l_OB.*m2;
t8 = t4+t5+t6+t7;
Corr_Joint_Sp = [-dth2.*t2.*t8.*(dth1.*2.0+dth2);dth1.^2.*t2.*t8;-dth4.*t3.*t8.*(dth3.*2.0+dth4);dth3.^2.*t3.*t8];
end
