function drOE2 = velocity_foot2(in1,in2)
%VELOCITY_FOOT2
%    drOE2 = VELOCITY_FOOT2(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    04-Nov-2024 14:09:58

dth3 = in1(7,:);
dth4 = in1(8,:);
l_AC = in2(17,:);
l_DE = in2(18,:);
l_OB = in2(16,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = cos(th3);
t3 = sin(th3);
t4 = th3+th4;
t5 = cos(t4);
t6 = sin(t4);
drOE2 = [dth3.*(l_AC.*t5+l_DE.*t2+l_OB.*t2)+dth4.*l_AC.*t5;dth3.*(l_AC.*t6+l_DE.*t3+l_OB.*t3)+dth4.*l_AC.*t6];
end
