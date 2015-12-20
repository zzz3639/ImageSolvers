function [ Mres, Index ] = Merge( pic, th, Izero)
%Modified in 2015.12.15 by ZHANG Haowen
%  Output:
%    Mres: Mres is the merged results, its length equals to number of merged groups.
%    Index: which merges to which, length equals to pic.
%  Input:
%    pic: [X,Y,I].
%    th: length within which 2 points are merged.
%    Izero: Points lower than this intensity are ignored.
n=size(pic,1);
Remove=pic(:,3)<Izero;
picR=zeros(n-sum(Remove),size(pic,2));
Index1=zeros(n,1);
u=(pic(:,3)>Izero);
Index1(find(u))=[1:length(find(u))];
picR=pic(u,:);

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
Index2=zeros(size(M,1),1);
Index2=M(v,4);
Index=zeros(size(Index1,1),1);
v=find(Index1);
Index(v,:)=Index2(Index1(v,:),:);

%re-arrange the index such that the point sequences before and after merge are consist.
[u,v1]=sort(Index);
[u,v2]=unique(u);
if Index(v1(v2(1)))==0
    S=v1(v2(2:end));
else
    S=v1(v2);
end
[u,vs]=sort(S);
IndexAns=zeros(size(Index1,1),1);
v=find(Index);
IndexAns(v,:)=vs(Index(v,:));
Index=IndexAns;
Mres=Mres(vs,:);

end

