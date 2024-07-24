load("parameter.mat");
w=1.4005;%入射波浪频率(s-1)
m_z_add=1335.535;%垂荡附加质量(kg)
I_x_add=6779.315;%纵摇附加转动惯量(kg·m2)
eta_z=656.3616;%垂荡兴波阻尼系数(N·s/m)
eta_x=151.4388;%纵摇兴波阻尼系数(N·m·s)
F_z=6250;%垂荡激励力振幅(N)
F_x=1230;%纵摇激励力矩振幅(N·m)
eta_d=10000;%直线阻尼系数(N·s/m)
fun=@(t,x)[x(3) x(4) (F_z*cos(w*t(1))+eta_d*abs(x(4)-x(3))^0.5*(x(4)-x(3))+...
    k_F*(x(2)-x(1))-eta_z*x(3)-pho*g*R_o^2*pi*x(1))/(m_o+m_z_add) ...
    -(eta_d*abs(x(4)-x(3))^0.5*(x(4)-x(3))+k_F*(x(2)-x(1)))/m_i]';
%x(1)为浮子位移,x(2)为振子位移,x(3)为浮子速度,x(4)为振子速度
t0=0;
x0=[0 0 0 0]';%初值
t_e=2*pi/w*40;%时长
step=5e-4;%步长
acc=10;%精度阶数
[X,T]=iAdams(fun,x0,t0,t_e,acc,step);%使用Adams内插法求解
T=T(1,:);
subplot(1,2,1)
plot(T,X(1,:),'r')
hold on
plot(T,X(1,:)-X(2,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('位移 m')
subplot(1,2,2)
plot(T,X(3,:),'r')
hold on
plot(T,X(3,:)-X(4,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('速度 m/s')
restep=0.2/step*(1:length(0:0.2:t_e)-1);
writematrix([T(1,restep)' X(1,restep)' X(2,restep)' X(3,restep)' X(4,restep)'],'result1-2.xlsx')
disp(X(:,[10 20 40 60 100]/step));