function q = q_twolegs(in1,in2)
%Q_TWOLEGS
%    Q = Q_TWOLEGS(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    11-Nov-2024 23:02:07

th1 = in1(3,:);
th2 = in1(4,:);
th3 = in1(5,:);
th4 = in1(6,:);
x = in1(1,:);
y = in1(2,:);
q = [x;y;th1;th2;th3;th4];
end