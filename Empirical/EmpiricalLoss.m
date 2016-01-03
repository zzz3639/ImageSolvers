function [Losses] = EmpiricalLoss( Img, Sigma, EstPara, PathPara, LambdaPath, EstimateLambda0 )
%EMPIRICALLOSS Summary of this function goes here
%Usage: [Losses] = EmpiricalLoss(Img,Sigma,[Lambda0,N_Est],[Nrepeats,MaxIteNode,IzeroNode,PzeroNode],LambdaPath, EstimateLambda0/nothing)
%  

%%% Initialize
Nrepeats=PathPara(1);
MaxIteNode=PathPara(2);
IzeroNode=PathPara(3);
PzeroNode=PathPara(4);
S1=size(Img,1);
S2=size(Img,2);
n=length(LambdaPath);
NumMol=floor(size(Img,1)*size(Img,2)/50/Sigma/Sigma*20);
bsize=4;
psfdecay=5;
LambdaForNoise=0.2;
MergeDist=0.1;
MolZero=0.5;
FoldRemain=0.99;
FilterFold=0.75;
MaxIte0=500;
MaxIteFast=5000;
MaxIteTunning=5000;
IzeroSolver=sum(sum(Img))*1e-6;
PzeroSolver=1e-5;


%%% Estimate molecule list and noise scale
[NoScale] = NoiseEst(Img, Sigma, EstPara(2), LambdaForNoise);

[Pic,NoB]=RunSolverTunning(Img, NumMol, EstPara(1), [Sigma,bsize,psfdecay], [MaxIte0,MaxIteFast,IzeroSolver,PzeroSolver], ...
[MaxIteTunning,IzeroSolver,PzeroSolver], [MergeDist*Sigma,MolZero], [MergeDist*Sigma,MolZero,FoldRemain]);

No=NoB/(S1+2*bsize)/(S2+2*bsize)*S1*S2;

%%% Generate new image series based on noise scale
[ImgSilicon] = SiliconImaging(Nrepeats, [S1,S2], Pic, No, Sigma,NoScale);

%%% Run ReverseLasso through LambdaPath and Calculate Empirical loss
Losses=zeros(n,1);
for t=1:Nrepeats
    [ResultsPath] = ReverseLasso(ImgSilicon{t}, NumMol, LambdaPath, [Sigma,bsize,psfdecay], ...
    [MaxIteNode,IzeroNode,PzeroNode], [MergeDist*Sigma,MolZero]);
    for i=1:n
        PicThis=ResultsPath{i}.pic;
        AlighThis=MolListAlign([PicThis(:,1:2),Ints2Prob(PicThis(:,3),FilterFold)],Pic);
        Recover=zeros(size(Pic,1),1);
        for j=1:size(AlighThis,1)
            Recover(AlighThis(j,7))=Recover(AlighThis(j,7))+AlighThis(j,3);
        end
        RecoverValue=0;
        for j=1:size(Pic,1)
            if Recover(j)>1
                RecoverValue=RecoverValue+1/Recover(j);
            else
                RecoverValue=RecoverValue+Recover(j);
            end
        end
        RecoverValue=RecoverValue/size(Pic,1);
        E=sum(AlighThis(:,3).*AlighThis(:,4))/sum(AlighThis(:,3));
        Losses(i)=Losses(i)+E/RecoverValue;
    end
end

Losses=Losses/Nrepeats;

if exist('EstimateLambda0')
    [u,v]=min(Losses);
    Lambda0=LambdaPath(v)+0.05;
    [Pic,NoB]=RunSolverTunning(Img, NumMol, EstPara(1), [Sigma,bsize,psfdecay], [MaxIte0,MaxIteFast,IzeroSolver,PzeroSolver], ...
    [MaxIteTunning,IzeroSolver,PzeroSolver], [MergeDist*Sigma,MolZero], [MergeDist*Sigma,MolZero,FoldRemain]);
    No=NoB/(S1+2*bsize)/(S2+2*bsize)*S1*S2;
    [ImgSilicon] = SiliconImaging(Nrepeats, [S1,S2], Pic, No, Sigma,NoScale);
    
    Losses=zeros(n,1);
    for t=1:Nrepeats
        [ResultsPath] = ReverseLasso(ImgSilicon{t}, NumMol, LambdaPath, [Sigma,bsize,psfdecay], ...
        [MaxIteNode,IzeroNode,PzeroNode], [MergeDist*Sigma,MolZero]);
        for i=1:n
            PicThis=ResultsPath{i}.pic;
            AlighThis=MolListAlign([PicThis(:,1:2),Ints2Prob(PicThis(:,3),FilterFold)],Pic);
            Recover=zeros(size(Pic,1),1);
            for j=1:size(AlighThis,1)
                Recover(AlighThis(j,7))=Recover(AlighThis(j,7))+AlighThis(j,3);
            end
            RecoverValue=0;
            for j=1:size(Pic,1)
                if Recover(j)>1
                    RecoverValue=RecoverValue+1/Recover(j);
                else
                    RecoverValue=RecoverValue+Recover(j);
                end
            end
            RecoverValue=RecoverValue/size(Pic,1);
            E=sum(AlighThis(:,3).*AlighThis(:,4))/sum(AlighThis(:,3));
            Losses(i)=Losses(i)+E/RecoverValue;
        end
    end
    Losses=Losses/Nrepeats;
end

end

