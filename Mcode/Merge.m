function [ Mres ] = Merge( pic, th, Izero)
%   Modified in 2015.06.20 by ZHANG Haowen
%   Detailed explanation goes here
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

end

