function keypoints = keypoints_leg(in1,in2)
%KEYPOINTS_LEG
%    KEYPOINTS = KEYPOINTS_LEG(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Nov-2024 13:27:31

l_AC = in2(17,:);
l_DE = in2(18,:);
l_OA = in2(15,:);
l_OB = in2(16,:);
th1 = in1(3,:);
th2 = in1(4,:);
x = in1(1,:);
y = in1(2,:);
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = l_OA.*t2;
t6 = l_OB.*t2;
t7 = cos(t4);
t8 = l_OA.*t3;
t9 = l_OB.*t3;
t10 = sin(t4);
t11 = l_AC.*t7;
t12 = l_AC.*t10;
t13 = -t5;
t14 = -t6;
t15 = -t11;
keypoints = reshape([x,y,t8+x,t13+y,t9+x,t14+y,t8+t12+x,t13+t15+y,t9+t12+x,t14+t15+y,t9+t12+x+l_DE.*t3,t14+t15+y-l_DE.*t2],[2,6]);
end