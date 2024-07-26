function P=p4(innum)
load("parameter.mat");
w=1.9806;%入射波浪频率(s-1)
m_z_add=1091.099;%垂荡附加质量(kg)
I_x_add=7142.493;%纵摇附加转动惯量(kg·m2)
eta_z=528.5018;%垂荡兴波阻尼系数(N·s/m)
eta_x=1655.909;%纵摇兴波阻尼系数(N·m·s)
F_z=1760;%垂荡激励力振幅(N)
F_x=2140;%纵摇激励力矩振幅(N·m)
eta_d=innum(1);%直线阻尼系数(N·s/m)
eta_c=innum(2);%曲线阻尼系数(N·m·s)
I_o=8289.43;%浮子的转动惯量(kg·m2)
%I_i=202.75+2433*(0.75+x(2)-x(1))^2;振子的转动惯量(kg·m2)
fun=@(t,x)[x(5) x(6) x(7) x(8) ...
    (F_z*cos(w*t(1))+eta_d*(x(6)-x(5))+k_F*(x(2)-x(1))...
    -eta_z*x(5)-pho*g*R_o^2*pi*x(1))/(m_o+m_z_add) ...
    -(eta_d*(x(6)-x(5))+k_F*(x(2)-x(1)))/m_i ...
    (F_x*cos(w*t(1))+eta_c*(x(8)-x(7))+k_M*(x(4)-x(3))...
    -eta_x*x(7)-M_p*x(3))/(I_o+I_x_add) ...
    -(eta_c*(x(8)-x(7))+k_M*(x(4)-x(3)))/(202.75+2433*(0.75+x(2)-x(1))^2)]';
%x(1)为浮子位移,x(2)为振子位移,x(3)为浮子角位移,x(4)为振子角位移
%x(5)为浮子速度,x(6)为振子速度,x(6)为浮子角速度,x(8)为振子角速度
t0=0;
x0=[0 0 0 0 0 0 0 0]';%初值
t_e=180;%时长
step=1e-4;%步长
acc=10;%精度阶数
[X,~]=iAdams(fun,x0,t0,t_e,acc,step);%使用Adams内插法求解
restep=int32((100:step:180)/step);
P=sum(abs(X(6,restep)-X(5,restep)).^2)*step*eta_d/80+...
    sum(abs(X(8,restep)-X(7,restep)).^2)*step*eta_c/80;