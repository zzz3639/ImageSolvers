function [ pic, no, mv, psfno] = FastSolver( Img, NumMol, Lambda, OptPara, Tolerance, ProcessPara )
%REVERSELASSO Summary of this function goes here
% Modified by Zhang Haowen in 2015.10.28
% Usage: [ pic, no, mv ] = FastSolver( Img, NumMol, Lambda, OptPara, Tolerance, ProcessPara )
% Input variables:
%   Img: m by n matrix, raw image
%   NumMol: number of initial molecules, before merge. Integer
%   Lambda: Lambda for the optimization.
%   OptPara: [Sigma,BoundarySize,PSFdecay] or [Sigma,BoundarySize,PSFdecay,NoiseGrid1,NoiseGrid2]
%      if OptPara=[Sigma,BoundarySize,PSFdecay] then background noise is treated as even.
%      Otherwise noise is modeled by first order b-spline base.
%   Tolerance: [MaxIte0, MaxIte, Izero, Pzero]
%      MaxIte0: Maximal iteration number of the first SparseGMM run
%      MaxIte: Maximal iteration number of the second SparseGMM run
%      Izero: minimal Intensity change.
%      Pzero: minimal Position change.
%   ProcessPara: [MergeDist, MolZero]
%      MergeDist: Molecules within MergeDist are treated as one
%      MolZero: Molecules below this intensity are ignored.
% Output variables:
%   Results:
%      Results under each lambda.
%   IniMv:
%      Initial molecules table before optimization.

%%% plant random seed
rng('shuffle');

%%% parameter initialize
n=NumMol;
sigma=OptPara(1);
bsize=OptPara(2);
bdecay=OptPara(3);
if length(OptPara)==3
    EvenNoise=1;
else
    EvenNoise=0;
    NoiseGridSize1=optPARA(4);
    NoiseGridSize2=optPARA(5);
end
MaxIte0=Tolerance(1);
MaxIte =Tolerance(2);
Izero=Tolerance(3);
Pzero=Tolerance(4);
MergeDist=ProcessPara(1);
MolZero=ProcessPara(2);

%%% define image size and reshape Img
bm=Img;
s1=size(bm,1);
s2=size(bm,2);
b=reshape(bm,s1*s2,1);
sbs=(s1+2*bsize)*(s2+2*bsize);

%%% define image positions list
cx1=[0:s1-1]';
cx1=repmat(cx1,1,s2);
cx2=[0:s2-1];
cx2=repmat(cx2,s1,1);
cxo=[reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];

%%% Initialize Molecule table
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
mvthis.pic=pic;
mvthis.no=no;

%%% Run Solver and get a primary result
if EvenNoise==1
    [picp,nop,mv]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte0,Izero,Pzero],bsize,bdecay,Lambda,[int64(10000)],mvthis);
    psfno=ones(s1+2*bsize,s2+2*bsize)/(s1+2*bsize)/(s2+2*bsize);
else
    [picp,nop,mv,psfno]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte0,Izero,Pzero],bsize,bdecay,[NoiseGridSize1,NoiseGridSize2],Lambda,[int64(10000)],mvthis);
end
picp=Merge(picp,MergeDist,MolZero);
mvthis.pic=picp;
mvthis.no=nop;
if EvenNoise==1
    [pic,no]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte,Izero,Pzero],bsize,bdecay,Lambda,[int64(10000)],mvthis);
else
    [pic,no]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte,Izero,Pzero],bsize,bdecay,[NoiseGridSize1,NoiseGridSize2],Lambda,[int64(10000)],mvthis);
end

end