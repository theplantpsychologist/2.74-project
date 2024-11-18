function J1 = jacobian_foot1(in1,in2)
%JACOBIAN_FOOT1
%    J1 = JACOBIAN_FOOT1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Nov-2024 23:02:07

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
t7 = l_AC.*t5;
t8 = l_AC.*t6;
J1 = reshape([1.0,0.0,0.0,1.0,t7+l_DE.*t2+l_OB.*t2,t8+l_DE.*t3+l_OB.*t3,t7,t8,0.0,0.0,0.0,0.0],[2,6]);
end