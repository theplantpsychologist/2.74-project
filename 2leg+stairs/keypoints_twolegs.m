function keypoints = keypoints_twolegs(in1,in2)
%KEYPOINTS_TWOLEGS
%    KEYPOINTS = KEYPOINTS_TWOLEGS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    21-Nov-2024 19:08:00

l_AC = in2(17,:);
l_DE = in2(18,:);
l_OA = in2(15,:);
l_OB = in2(16,:);
th1 = in1(1,:);
th2 = in1(2,:);
th3 = in1(3,:);
th4 = in1(4,:);
t2 = cos(th1);
t3 = cos(th3);
t4 = sin(th1);
t5 = sin(th3);
t6 = th1+th2;
t7 = th3+th4;
t8 = l_OA.*t2;
t9 = l_OB.*t2;
t10 = l_OA.*t3;
t11 = l_OB.*t3;
t12 = cos(t6);
t13 = cos(t7);
t14 = l_OA.*t4;
t15 = l_OB.*t4;
t16 = l_OA.*t5;
t17 = l_OB.*t5;
t18 = sin(t6);
t19 = sin(t7);
t20 = l_AC.*t12;
t21 = l_AC.*t13;
t22 = l_AC.*t18;
t23 = l_AC.*t19;
t24 = -t8;
t25 = -t9;
t26 = -t10;
t27 = -t11;
t28 = -t20;
t29 = -t21;
keypoints = reshape([t14,t24,t15,t25,t14+t22,t24+t28,t15+t22,t25+t28,t15+t22+l_DE.*t4,t25+t28-l_DE.*t2,t16,t26,t17,t27,t16+t23,t26+t29,t17+t23,t27+t29,t17+t23+l_DE.*t5,t27+t29-l_DE.*t3],[2,10]);
end