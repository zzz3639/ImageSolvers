function [Losses] = EmpiricalLoss( Img, Sigma, EstPara, PathPara, LambdaPath )
%EMPIRICALLOSS Summary of this function goes here
%Usage: [Losses] = EmpiricalLoss(Img,Sigma,[Lambda0,N_Est],[Nrepeats,MaxIteNode,IzeroNode,PzeroNode],LambdaPath)

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
bdecay=5;
MergeDist=0.1;
MolZero=0.5;
FilterFold=0.75;

%%% Estimate molecule list and noise scale
[NoScale, Pic, No] = NoiseEst(Img, Sigma, EstPara(2), EstPara(1));

%%% Generate new image series based on noise scale
[ImgSilicon] = SiliconImaging(Nrepeats, [S1,S2], Pic, No, Sigma,NoScale);

%%% Run ReverseLasso through LambdaPath and Calculate Empirical loss
Losses=zeros(n,1);
for t=1:Nrepeats
    [ResultsPath] = ReverseLasso(ImgSilicon{t}, NumMol, LambdaPath, [Sigma,bsize,bdecay], ...
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

