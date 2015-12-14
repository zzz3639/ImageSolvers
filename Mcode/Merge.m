function [ Mres, Index ] = Merge( pic, th, Izero)
% Modified in 2015.12.04 by ZHANG Haowen
%  Output:
%    Mres: Mres is the merged results, its length equals to number of merged groups.
%    Index: which merges to which, length equals to pic.
n=size(pic,1);
Remove=pic(:,3)<Izero;
picR=zeros(n-sum(Remove),size(pic,2));
j=1;
for i=1:n
    if pic(i,3)<Izero;
        continue;
    end
    picR(j,:)=pic(i,:);
    j=j+1;
end

n=size(picR,1);
pic=picR;

M=MergeMex(pic,th);
M(:,3:4)=M(:,3:4)+1;

u=sort(M(:,4));
u=unique(u);
nr=length(u);

Mres=zeros(nr,3);
for i=1:n
    Mres(M(i,4),1:2)=Mres(M(i,4),1:2)+pic(M(i,3),1:2)*pic(M(i,3),3);
    Mres(M(i,4),3)=Mres(M(i,4),3)+pic(M(i,3),3);
end

Mres(:,1:2)=Mres(:,1:2)./repmat(Mres(:,3),1,2);
[u,v]=sort(M(:,3));
Index=zeros(size(M,1),1);
Index=M(v,4);

end

