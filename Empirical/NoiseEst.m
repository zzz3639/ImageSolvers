function [NoScale, Pic, No] = NoiseEst(Img, Sigma, NTrial, Lambda0)
%NOISEEST Summary of this function goes here
%   Detailed explanation goes here

%%% Define some useful constants  
MaxIte0=500;
MaxIteFast=5000;
MaxIteTunning=5000;
IzeroSolver=sum(sum(Img))*1e-5;
PzeroSolver=1e-4;
IzeroTrial=sum(sum(Img))*1e-5;
PzeroTrial=1e-5;
bsize=4;
bdecay=5;
MergeDist=0.1;
MolZero=0.5;
FoldRemain=0.9999;

%%% Initialize
rng('shuffle');
S1=size(Img,1);
S2=size(Img,2);
NumMol=floor(size(Img,1)*size(Img,2)/50/Sigma/Sigma*20);
SearchMin=0;
SearchMax=0.5;
SearchBar=0.02;

%%% Run RunSolverTunning.m and get squared error
[Pic,No,Error]=ErrorEst(Img, NumMol, Lambda0, [Sigma,bsize,bdecay], ...
[MaxIte0,MaxIteFast,IzeroSolver,PzeroSolver], [MaxIteTunning,IzeroSolver,PzeroSolver], ...
[MergeDist*Sigma,MolZero], [MergeDist*Sigma,MolZero,FoldRemain]);

%%% Estimate NoiseScale by binary search and SiliconImaging
while SearchMax-SearchMin>SearchBar
    ObNoise=(SearchMax+SearchMin)/2;
    [ImgSilicon] = SiliconImaging(NTrial, [S1,S2], Pic, No, Sigma,ObNoise);
    ErrorThisAve=0;
    for i=1:NTrial
        [PicThis,NoThis,ErrorThis]=ErrorEst(ImgSilicon{i},NumMol, Lambda0, [Sigma,bsize,bdecay], ...
        [MaxIte0,MaxIteFast,IzeroTrial,PzeroTrial], [MaxIteTunning,IzeroTrial,PzeroTrial], ...
        [MergeDist*Sigma,MolZero], [MergeDist*Sigma,MolZero,FoldRemain]);
        ErrorThisAve=ErrorThisAve+ErrorThis;
    end
    ErrorThisAve=ErrorThisAve/NTrial;
    if Error>ErrorThisAve
        SearchMin=ObNoise;
    else
        SearchMax=ObNoise;
    end
end

NoScale=(SearchMax+SearchMin)/2;

end

function [Pic,No,Error]=ErrorEst(Img,NumMol,Lambda0,OptPara,StopPara1,StopPara2,ProcessingPara1,ProcessingPara2)
    S1=size(Img,1);
    S2=size(Img,2);
    bsize=OptPara(2);
    Sigma=OptPara(1);
    cx1=[0:S1-1]';
    cx1=repmat(cx1,1,S2);
    cx2=[0:S2-1];
    cx2=repmat(cx2,S1,1);
    cx=[reshape(cx2,S1*S2,1),reshape(cx1,S1*S2,1)];
    [Pic,NoB]=RunSolverTunning(Img,NumMol,Lambda0,OptPara,StopPara1,StopPara2,ProcessingPara1,ProcessingPara2);

    n=size(Pic,1);
    No=NoB/(S1+2*bsize)/(S2+2*bsize)*S1*S2;
    ImgReal=ones(S1*S2,1)*No/S1/S2;
    for i=1:n
        ImgReal=ImgReal+Pic(i,3)*PSF(Pic(i,1:2),cx,Sigma);
    end
    ImgReal=reshape(ImgReal,S1,S2);
    Error=sum(sum((Img-ImgReal).^2));
end

function b2=Noise(s1,s2,PhotonNum)
    b2=mnrnd(PhotonNum,repmat(1/s1/s2,s1*s2,1));
end

%PSF function
function V=PSF(mu,cx,sig)
    V=normpdf(cx(:,1),mu(1),sig).*normpdf(cx(:,2),mu(2),sig);
    V=V/sum(V);
end


