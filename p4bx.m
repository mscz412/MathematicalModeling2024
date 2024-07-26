P=zeros(25,25);
parfor i=1:25
    for j=1:25
        P(i,j)=p4([4000*i 4000*j]);
    end
end
mesh(4000*(1:25),4000*(1:25),P)
xlabel('直线阻尼系数')
ylabel('旋转阻尼系数')
zlabel('P(w)')