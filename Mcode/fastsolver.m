function [ pic, no, mv ] = fastsolver( Img, NumMol, Lambda, OptPara, Tolerance, ProcessPara )
%REVERSELASSO Summary of this function goes here
% Modified by Zhang Haowen in 2015.08.19
% Input variables:
%   Img: m by n matrix, raw image
%   NumMol: number of initial molecules, before merge. Integer
%   Lambda: Lambda for the optimization.
%   OptPara: [Sigma, BoundarySize, PSFdecay];
%   Tolerance: [MaxIte0, MaxIte, Izero, Pzero]. 
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
no=[0.5];
scale=sum(sum(b))+(sbs-s1*s2)*median(b);
no=no*scale;
pic(:,3)=pic(:,3)*scale;
mvthis.pic=pic;
mvthis.no=no;

%%% Run Solver and get a primary result
[picp,nop,mv]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte0,Izero,Pzero],bsize,bdecay,Lambda,[int64(10000)],mvthis);
picp=Merge(picp,MergeDist,MolZero);
mvthis.pic=picp;
mvthis.no=nop;
[pic,no]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte,Izero,Pzero],bsize,bdecay,Lambda,[int64(10000)],mvthis);

end

