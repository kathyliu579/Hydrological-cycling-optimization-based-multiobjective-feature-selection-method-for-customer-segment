function [pEva] = createPossibility(pop,nPop)
dis=zeros(1,nPop);
d=dis;
pEva=d;
for i=1:nPop
    for j=1:nPop
        dis(j)=sum(abs(pop(i).Position-pop(j).Position));
    end
    d(i)=sum(dis);
end

averageDis=sum(d)/nPop;
for i=1:nPop
    pEva(i)=exp(-(d(i)/averageDis)^3);
end
end



