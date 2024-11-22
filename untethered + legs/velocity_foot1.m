function drOE1 = velocity_foot1(in1,in2)
%VELOCITY_FOOT1
%    drOE1 = VELOCITY_FOOT1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    21-Nov-2024 19:07:38

dth1 = in1(9,:);
dth2 = in1(10,:);
dx = in1(7,:);
dy = in1(8,:);
l_AC = in2(17,:);
l_DE = in2(18,:);
l_OB = in2(16,:);
th1 = in1(3,:);
th2 = in1(4,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = cos(t4);
t6 = sin(t4);
drOE1 = [dx+dth1.*(l_AC.*t5+l_DE.*t2+l_OB.*t2)+dth2.*l_AC.*t5;dy+dth1.*(l_AC.*t6+l_DE.*t3+l_OB.*t3)+dth2.*l_AC.*t6];
end
