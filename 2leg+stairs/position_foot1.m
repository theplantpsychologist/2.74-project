function rOE1 = position_foot1(in1,in2)
%POSITION_FOOT1
%    rOE1 = POSITION_FOOT1(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    21-Nov-2024 19:07:59

l_AC = in2(17,:);
l_DE = in2(18,:);
l_OB = in2(16,:);
th1 = in1(1,:);
th2 = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
rOE1 = [l_DE.*t3+l_OB.*t3+l_AC.*sin(t4);-l_DE.*t2-l_OB.*t2-l_AC.*cos(t4)];
end
