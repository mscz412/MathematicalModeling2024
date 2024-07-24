function P=p2_1(r)
load("parameter.mat");
w=2.2143;%入射波浪频率(s-1)
m_z_add=1165.992;%垂荡附加质量(kg)
I_x_add=7131.29;%纵摇附加转动惯量(kg·m2)
eta_z=167.8395;%垂荡兴波阻尼系数(N·s/m)
eta_x=2992.724;%纵摇兴波阻尼系数(N·m·s)
F_z=4890;%垂荡激励力振幅(N)
F_x=2560;%纵摇激励力矩振幅(N·m)
eta_d=r;%直线阻尼系数(N·s/m)
fun=@(t,x)[x(3) x(4) (F_z*cos(w*t(1))+eta_d*(x(4)-x(3))+...
    k_F*(x(2)-x(1))-eta_z*x(3)-pho*g*R_o^2*pi*x(1))/(m_o+m_z_add) ...
    -(eta_d*(x(4)-x(3))+k_F*(x(2)-x(1)))/m_i]';
%x(1)为浮子位移,x(2)为振子位移,x(3)为浮子速度,x(4)为振子速度
t0=0;
x0=[0 0 0 0]';%初值
t_e=180;%时长
step=1e-4;%步长
acc=10;%精度阶数
[X,T]=iAdams(fun,x0,t0,t_e,acc,step);%使用Adams内插法求解
restep=int32((100:step:180)/step);
P=sum(abs(X(4,restep)-X(3,restep)).^2)*step*eta_d/80;