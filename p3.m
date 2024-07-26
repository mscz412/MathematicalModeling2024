load("parameter.mat");
w=1.7152;%入射波浪频率(s-1)
m_z_add=1028.876;%垂荡附加质量(kg)
I_x_add=7001.914;%纵摇附加转动惯量(kg·m2)
eta_z=683.4558;%垂荡兴波阻尼系数(N·s/m)
eta_x=654.3383;%纵摇兴波阻尼系数(N·m·s)
F_z=3640;%垂荡激励力振幅(N)
F_x=1690;%纵摇激励力矩振幅(N·m)
eta_d=10000;%直线阻尼系数(N·s/m)
eta_c=1000;%曲线阻尼系数(N·m·s)
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
t_e=2*pi/w*40;%时长
step=1e-4;%步长
acc=10;%精度阶数
[X,T]=iAdams(fun,x0,t0,t_e,acc,step);%使用Adams内插法求解
T=T(1,:);
subplot(2,2,1)
plot(T,X(1,:),'r')
hold on
plot(T,X(1,:)-X(2,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('位移 m')
subplot(2,2,2)
plot(T,X(3,:),'r')
hold on
plot(T,X(3,:)-X(4,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('角位移 rad')
subplot(2,2,3)
plot(T,X(5,:),'r')
hold on
plot(T,X(5,:)-X(6,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('速度 m/s')
subplot(2,2,4)
plot(T,X(7,:),'r')
hold on
plot(T,X(7,:)-X(8,:),'b')
legend('浮子','振子')
xlabel('时间 s')
ylabel('角速度 rad/s')
restep=int32(0.2/step*(1:length(0:0.2:t_e)-1));
writematrix([T(1,restep)' X(1,restep)' X(2,restep)' X(3,restep)' X(4,restep)' ...
    X(5,restep)' X(6,restep)' X(7,restep)' X(8,restep)'],'result3.xlsx')
writematrix([[10 20 40 60 100];X(:,[10 20 40 60 100]/step)],'p3.xlsx')