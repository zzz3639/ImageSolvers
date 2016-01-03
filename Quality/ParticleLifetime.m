function [ Lower, Upper ] = ParticleLifetime( LinkPaths, LambdaPaths )
%Modified by ZHANG Haowen in 2015.12.17
%   Detailed explanation goes here

n=length(LambdaPaths);
Lower=cell(n,1);
upper=cell(n,1);
FullIndex=cell(n-1,1);
FullIndexReverse=cell(n-1,1);

% find non-merge non-decay links
for i=1:n-1
    L=LinkPaths{i+1};
    v1=find(L);
    [SL,v2]=sort(L(v1));
    DSL=([1;SL(2:end)-SL(1:end-1)]~=0);
    ZDSL=find(DSL==0);
    DSL(ZDSL-1)=0;
    v3=find(DSL==0);
    L(v1(v2(v3)))=0;
    FullIndex{i}=L;
    Nm=max(LinkPaths{i+1});
    Lr=zeros(Nm,1);
    v=find(L);
    Lr(L(v))=v;
    FullIndexReverse{i}=Lr;
end

% calculate Lower
Lower{1}=zeros(length(LinkPaths{2}),1);
Live=ones(length(LinkPaths{2}),1);
for i=2:n
    L=FullIndex{i-1};
    Lr=FullIndexReverse{i-1};
    LiveNew=i*ones(length(Lr),1);
    v=find(Lr);
    LiveNew(v)=Live(Lr(v));
    LowerThis=zeros(length(Lr),1);
    LowerThis(v)=i-Live(Lr(v));
    Lower{i}=LowerThis;
    Live=LiveNew;
end

% calculate Upper
Upper{n}=zeros(length(FullIndexReverse{n-1}),1);
Live=n*ones(length(FullIndexReverse{n-1}),1);
for i=n-1:-1:1
    L=FullIndex{i};
    Lr=FullIndexReverse{i};
    LiveNew=i*ones(length(L),1);
    v=find(L);
    LiveNew(v)=Live(L(v));
    UpperThis=zeros(length(L),1);
    UpperThis(v)=Live(L(v))-i;
    Upper{i}=UpperThis;
    Live=LiveNew;
end


end

