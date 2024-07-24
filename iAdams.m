function [Y,t]=iAdams(f,y0,t0,b,acc,step)
    sizey=length(y0);
    stept=step*ones(sizey,1);
    Y=y0;
    t=t0*ones(sizey,1);
    fs=zeros(sizey,acc);
    fs(:,1)=f(t,y0);
    y=y0;
    for i=1:acc-1
        y(:,i+1)=y(:,i)-step*f(t-stept*(i-1),y(:,i));
        last=y(:,i+1);
        while last~=y(:,i+1)
            last=y(:,i+1);
            y(:,i+1)=y(:,i)-step*f(t-stept*i,y(:,i+1));
        end
        fs(:,i+1)=f(t-stept*i,y(:,i+1));
    end
    i=1;
    B=getb(acc+1);
    while b-t(1,i)>=step
        t(:,i+1)=t(:,i)+stept;
        Y(:,i+1)=Y(:,i)+step*(B(acc,1:acc)*fs(:,1:acc)')';
        fs=[f(t(:,i),Y(:,i)) fs];
        last=Y(i+1)+1;
        while last~=Y(:,i+1)
            last=Y(:,i+1);
            Y(:,i+1)=Y(:,i)+step*(B(acc+1,1:acc+1)*fs(:,1:acc+1)')';
            fs(:,1)=f(t(:,i+1),Y(:,i+1));
        end
        i=i+1;
    end
    if b-t(1,i)>=1e-10
        Y(:,i+1)=Y(:,i)+(b-t(1,i))*(B(acc,1:acc)*fs(:,1:acc)')';
        t(:,i+1)=b*ones(sizey,1);
    end
end
function A=geta(acc)
    A=1;
    K=1/2;
    for i=1:acc-1
        A=[A -A*K'];
        K=[1/(i+2) K];
    end
end
function B=getb(acc)
    A=geta(acc+1);
    B=zeros(acc,acc);
    for i=1:acc
        for k=i:acc
            C=ones(1,k-i+1);
            for j=i:k+1
                C(j-i+1)=factorial(j-1)/factorial(i-1)/factorial(j-i);
            end
            B(k,i)=A(i:k+1)*C';
        end
    end
    for i=1:acc
        B(:,i)=(2*mod(i,2)-1)*B(:,i);
    end
end  