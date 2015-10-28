function [ pic,no,mv,psfno ] = RunSolver( bm, n, sigma, StopPara, optPARA, lambda, SavePara, mv0 )
%  Modified in 2015.10.28, do the same thing as EMSparseSmooth
%    [pic,no,mv,Time] = EMboundary(bm, n, sigma, [maxite,Izero,Pzero], [bsize, bdecay] or [bsize,bdecay,NoiseGridSize1,NoiseGridSize2], lambda, save/kframes, mv0 or nothing);  
%    If OptPARA=[bsize, bdecay], then this code treat the background noise as even noise.
%    Otherwise noise is modeled by first order b-spline base.

%%% Maximal iteration number
iteration=StopPara(1);
%%% Stop when likelihood decrease < Lzero
Lzero=StopPara(2);

%%% plant random seed
rng('shuffle');

%%% define image size and reshape bm
bsize=optPARA(1);
bdecay=optPARA(2);
if length(optPARA)==2
    EvenNoise=1;
else
    EvenNoise=0;
    NoiseGridSize1=optPARA(3);
    NoiseGridSize2=optPARA(4);
end
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
    if EvenNoise==1
        no=[0.5];
    else
        numg1=floor((s1-1-floor(NoiseGridSize1/2))/NoiseGridSize1)+1;
        numg2=floor((s2-1-floor(NoiseGridSize2/2))/NoiseGridSize2)+1;
        NumNo=numg1*numg2;
        no=0.5/NumNo*ones(NumNo,1);
    end
    scale=sum(sum(b))+(sbs-s1*s2)*median(b);
    no=no*scale;
    pic(:,3)=pic(:,3)*scale;
    mv0.pic=pic;
    mv0.no=no;
end

if EvenNoise==1
    [pic,no,mv]=EMSparseSmoothMex(bm,n,sigma,StopPara,bsize,bdecay,lambda,[int64(SavePara(1))],mv0);
    psfno=ones(s1+2*bsize,s2+2*bsize)/(s1+2*bsize)/(s2+2*bsize);
else
    [pic,no,mv,psfno]=EMSparseUnevenMex(bm,n,sigma,StopPara,bsize,bdecay,[NoiseGridSize1,NoiseGridSize2],lambda,[int64(SavePara(1))],mv0);
end
end

