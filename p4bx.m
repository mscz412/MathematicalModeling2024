P=zeros(25,25);
parfor i=1:25
    for j=1:25
        P(i,j)=p2_2([4000*i 4000*j]);
    end
end