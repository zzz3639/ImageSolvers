function [ NoScale ] = NoiseEst(Img, Sigma, NTrial, Lambda0)
%NOISEEST Summary of this function goes here
%   Detailed explanation goes here

%%% Define some useful constants  
MaxIte0=500;
MaxIteFast=5000;
MaxIteTunning=5000;
Izero=sum(sum(Img))*1e-5;
Pzero=sum(sum(Img))*1e-4;
bsize=4;
bdecay=5;
MergeDist=0.1;
MolZero=0.5;
FoldRemain=0.9999;

%%% Initialize
rng('shuffle');
S1=size(Img,1);
S2=size(Img,2);
NumMol=floor(size(Img,1)*size(Img,2)/50*20);
cx1=[0:S1-1]';
cx1=repmat(cx1,1,S2);
cx2=[0:S2-1];
cx2=repmat(cx2,S1,1);
cx=[reshape(cx2,S1*S2,1),reshape(cx1,S1*S2,1)];
SearchMin=0;
SearchMax=0.5;
SearchBar=0.02;

%%% Run RunSolverTunning.m and get squared error
[Pic,NoB]=RunSolverTunning(Img, NumMol, Lambda0, [Sigma,bsize,bdecay], ...
[MaxIte0,MaxIteFast,Izero,Pzero], [MaxIteTunning,Izero,Pzero], ...
[MergeDist*Sigma,MolZero], [MergeDist*Sigma,MolZero,FoldRemain]);

n=size(Pic,1);
No=NoB/(S1+2*bsize)/(S2+2*bsize)*S1*S2;
ImgReal=ones(S1*S2,1)*No/S1/S2;
for i=1:n
    ImgReal=ImgReal+PSF(Pic(i,1:2),cx,Sigma);
end
ImgReal=reshape(ImgReal,S1,S2);
Error2=sum(sum((Img-ImgReal).^2));

%%% Estimate NoiseScale by binary search and SiliconImaging
while SearchMax-SearchMin>SearchBar
    ObNoise=(SearchMax+SearchMin)/2;
    [ImgSilicon] = SiliconImaging(NTrial, [S1,S2], Pic, No, Sigma);
    for i=1:NTrial
    end
end

end

function b2=Noise(s1,s2,PhotonNum)
    b2=mnrnd(PhotonNum,repmat(1/s1/s2,s1*s2,1));
end

%PSF function
function V=PSF(mu,cx,sig)
    V=normpdf(cx(:,1),mu(1),sig).*normpdf(cx(:,2),mu(2),sig);
    V=V/sum(V);
end


