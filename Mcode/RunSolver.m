function [ pic,no,mv ] = RunSolver( bm, n, sigma, StopPara, optPARA, lambda, SavePara, mv0 )
%    Modified in 2015.06.14, do the same thing as EMSparseSmooth
%    [pic,no,mv,Time] = EMboundary(bm, n, sigma, [maxite,Izero,Pzero], [bsize, bdecay], lambda, save/kframes, mv0/nothing);  

%%% Maximal iteration number
iteration=StopPara(1);
%%% Stop when likelihood decrease < Lzero
Lzero=StopPara(2);

%%% plant random seed
rng('shuffle');

%%% define image size and reshape bm
bsize=optPARA(1);
bdecay=optPARA(2);
s1=size(bm,1);
s2=size(bm,2);
b=reshape(bm,s1*s2,1);

%%% constants
sbs=(s1+2*bsize)*(s2+2*bsize);

%%% define image positions list
cx1=[0:s1-1]';
cx1=repmat(cx1,1,s2);
cx2=[0:s2-1];
cx2=repmat(cx2,s1,1);
cxo=[reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];

%%% Initialize Molecule table
if exist('mv0')
else
    Dtemp=mnrnd(n,b'/sum(b));
    Mtemp=zeros(n,1);
    [v]=find(Dtemp);
    k=1;
    for i=v
        Mtemp(k:k-1+Dtemp(i),:)=i;
        k=k+Dtemp(i);
    end
    pic=[cxo(Mtemp,1)+rand(n,1)-0.5,cxo(Mtemp,2)+rand(n,1)-0.5,ones(n,1)/(2*n)];
    no=[0.5];
    scale=sum(sum(b))+(sbs-s1*s2)*median(b);
    no=no*scale;
    pic(:,3)=pic(:,3)*scale;
    mv0.pic=pic;
    mv0.no=no;
end

[pic,no,mv]=EMSparseSmoothMex(bm,n,sigma,StopPara,bsize,bdecay,lambda,[int64(SavePara(1))],mv0);

end

