function A = A_leg(in1,in2)
%A_leg
%    A = A_leg(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    20-Nov-2024 15:01:32

I1 = in2(5,:);
I2 = in2(6,:);
I3 = in2(7,:);
I4 = in2(8,:);
Ir = in2(9,:);
N = in2(10,:);
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
t2 = cos(th1);
t3 = sin(th1);
t4 = th1+th2;
t5 = N.^2;
t6 = l_AC.^2;
t7 = l_A_m3.^2;
t8 = l_B_m2.^2;
t9 = l_O_m1.^2;
t10 = m2.*4.0;
t11 = m3.*4.0;
t12 = m4.*4.0;
t13 = Ir.*N;
t14 = l_C_m4.*t2;
t15 = l_OA.*t2;
t16 = l_OB.*t2;
t17 = cos(t4);
t18 = l_C_m4.*t3;
t19 = l_OA.*t3;
t20 = l_OB.*t3;
t21 = sin(t4);
t22 = l_O_m1.*m1.*t2;
t23 = l_O_m1.*m1.*t3;
t24 = Ir.*t5;
t51 = m1+t10+t11+t12;
t25 = t14.*4.0;
t26 = t15.*4.0;
t27 = t16.*4.0;
t28 = t18.*4.0;
t29 = t19.*4.0;
t30 = t20.*4.0;
t31 = t17.^2;
t32 = t21.^2;
t33 = l_AC.*t17;
t34 = l_A_m3.*t17;
t35 = l_B_m2.*t17;
t36 = l_AC.*t21;
t37 = l_A_m3.*t21;
t38 = l_B_m2.*t21;
t39 = t33.*4.0;
t40 = t34.*4.0;
t41 = t35.*4.0;
t42 = t36.*4.0;
t43 = t37.*4.0;
t44 = t38.*4.0;
t45 = m4.*t33.*2.0;
t46 = m3.*t34.*2.0;
t47 = m2.*t35.*2.0;
t48 = m4.*t36.*2.0;
t49 = m3.*t37.*2.0;
t50 = m2.*t38.*2.0;
t52 = t15+t34;
t53 = t16+t35;
t54 = t19+t37;
t55 = t20+t38;
t60 = t14+t15+t33;
t61 = t18+t19+t36;
t56 = t26+t40;
t57 = t27+t41;
t58 = t29+t43;
t59 = t30+t44;
t62 = t25+t26+t39;
t63 = t37.*t54.*2.0;
t64 = t38.*t55.*2.0;
t67 = t28+t29+t42;
t70 = t34.*t52.*2.0;
t71 = t35.*t53.*2.0;
t72 = t33.*t60.*2.0;
t73 = t36.*t61.*2.0;
t76 = t45+t46+t47;
t77 = t48+t49+t50;
t65 = (m3.*t56)./2.0;
t66 = (m2.*t57)./2.0;
t68 = (m3.*t58)./2.0;
t69 = (m2.*t59)./2.0;
t74 = (m4.*t67)./2.0;
t75 = (m4.*t62)./2.0;
t78 = t63+t70;
t79 = t64+t71;
t82 = t72+t73;
t80 = (m3.*t78)./2.0;
t81 = (m2.*t79)./2.0;
t83 = (m4.*t82)./2.0;
t84 = t23+t68+t69+t74;
t85 = t22+t65+t66+t75;
t86 = I2+I3+t13+t80+t81+t83;
A = reshape([t51,0.0,t85,t76,0.0,t51,t84,t77,t85,t84,I1+I2+I3+I4+Ir+t24+(m1.*(t2.^2.*t9.*2.0+t3.^2.*t9.*2.0))./2.0+(m3.*(t52.^2.*2.0+t54.^2.*2.0))./2.0+(m2.*(t53.^2.*2.0+t55.^2.*2.0))./2.0+(m4.*(t60.^2.*2.0+t61.^2.*2.0))./2.0,t86,t76,t77,t86,I2+I3+t24+(m4.*(t6.*t31.*2.0+t6.*t32.*2.0))./2.0+(m3.*(t7.*t31.*2.0+t7.*t32.*2.0))./2.0+(m2.*(t8.*t31.*2.0+t8.*t32.*2.0))./2.0],[4,4]);
end
