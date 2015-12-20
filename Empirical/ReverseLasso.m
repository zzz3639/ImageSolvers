function [ Results, IniMv, LinkPath ] = ReverseLasso( Img, NumMol, LambdaPath, OptPara, Tolerance, ProcessPara )
%REVERSELASSO Summary of this function goes here
% Modified by Zhang Haowen in 2015.07.10
% Input variables:
%   Img: m by n matrix, raw image
%   NumMol: number of initial molecules, before merge.
%   LambdaPath: Lambdas we want to try. Example: [0.02:0.02:0.3];
%   OptPara: [Sigma, BoundarySize, PSFdecay];
%   Tolerance: [MaxIteNode, Izero, Pzero]. 
%      MaxIteNode: Maximal iteration number of each lambda value in LambdaPath
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
%   LinkPath:
%      Record merge events through the LambdaPath.

%%% plant random seed
rng('shuffle');

%%% parameter initialize
n=NumMol;
sigma=OptPara(1);
bsize=OptPara(2);
bdecay=OptPara(3);
MaxIte=Tolerance(1);
Izero=Tolerance(2);
Pzero=Tolerance(3);
MergeDist=ProcessPara(1);
MolZero=ProcessPara(2);

%%% define image size and reshape Img
bm=Img;
s1=size(bm,1);
s2=size(bm,2);
b=reshape(bm,s1*s2,1);
sbs=(s1+2*bsize)*(s2+2*bsize)

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

%%% Define output variable
Results=cell(length(LambdaPath),1);
LinkPath=cell(length(LambdaPath),1);

%%% Run Sparse Solver with each value in LambdaPath
mvthis.pic=pic;
mvthis.no=no;
IniMv=mvthis;
for l=1:length(LambdaPath)
    lambda=LambdaPath(l);
    [picthis0,nothis]=EMSparseSmoothMex(bm,size(mvthis.pic,1),sigma,[MaxIte,Izero,Pzero],bsize,bdecay,lambda,[int64(10000)],mvthis);
    [picthis,indexthis]=Merge(picthis0,MergeDist,MolZero);
    mvthis.pic=picthis;
    mvthis.no=nothis;
    Results{l}=mvthis;
    LinkPath{l}=indexthis;
end

end

