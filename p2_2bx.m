P=zeros(25,25);
parfor i=1:25
    for j=1:25
        P(i,j)=p2_2([4000*i 0.04*j]);
    end
end
xlabel('直线阻尼系数')
ylabel('幂指数')
zlabel('P(w)')