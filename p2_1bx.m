P=zeros(1,100);
parfor i=1:100
    P(i)=p2_1(1000*i);
end
plot(1000*(1:100))
xlabel('直线阻尼系数')
ylabel('P(w)')