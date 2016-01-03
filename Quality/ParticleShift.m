function [ shift ] = ParticleShift( Rlasso, LinkPaths )
% Usage: [ shift ] = ParticleShift( Rlasso, LinkPaths )
%   Detailed explanation goes here

N=length(Rlasso);
shift=cell(N,1);
for i=1:N
    shift{i}=zeros(size(Rlasso{i}.pic,1),2);
end
for i=1:N-1
    PicF=Rlasso{i}.pic;
    PicT=Rlasso{i+1}.pic;
    L=LinkPaths{i+1};
    [u,v]=sort(L);
    S=[u;u(end)+1]-[-1;u];
    s=find(S);
    if L(v(s(1)))==0
        j=2;
    else
        j=1;
    end
    while j<length(s)
        f=s(j);
        t=s(j+1)-1;
        PT=PicT(u(t),:);
        PF=PicF(v(f:t),:);
        DF=sqrt((PF(:,1)-PT(1,1)).^2+(PF(:,2)-PT(1,2)).^2);
        shift{i}(v(f:t),2)=DF;
        shift{i+1}(u(t),1)=(DF'*PF(:,3))/sum(PF(:,3));
        j=j+1;
    end
    if L(v(s(1)))==0
        f=s(1);
        t=s(2)-1;
        shift{i}(v(f:t),2)=max(shift{i}(:,2));
    end
end

end

