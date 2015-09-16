function [ Pic ] = RunSeries( Sery, sig, Lambda0, NumTryNoise, NumTryLasso, LambdaPath, Skip )
%RUNSERIES Summary of this function goes here
%   Detailed explanation goes here
N=size(Sery,3);
MaxIteNode=2000;
Izero=1e-5;
Pzero=1e-4;
MaxIte0=5000;
Izero0=1e-6;
Pzero0=1e-5;


Pic=cell(N,1);
OptL=zeros(ceil(N/Skip),1);
k=1;
for t=1:N
    if mod(t-1,Skip)~=0
        continue;
    end
    Img=Sery(:,:,t);
    [Losses] = EmpiricalLoss(Sery(:,:,t),sig,[Lambda0,NumTryNoise],[NumTryLasso,MaxIteNode,Izero*sum(sum(Img)),Pzero],LambdaPath);
    [u,v]=min(Losses);
    OptL(k)=LambdaPath(v);
    k=k+1;
end

Norm=[0;ones(size(OptL,1)-1,1)]+ones(size(OptL))+[ones(size(OptL,1)-1,1);0];
OptL=[OptL(2:end,:);0]+OptL+[0;OptL(1:end-1,:)];
OptL=OptL./Norm;

for t=1:N
    Img=Sery(:,:,t);
    L=OptL(ceil(t/Skip));
    [pic,no]=RunSolver( Img, floor(size(Img,1)*size(Img,2)/7/7/sig/sig*20), sig, [MaxIte0,Izero0*sum(sum(Img)),Pzero0], [4,5], L, int64(10000) );
    pic=PostRun(pic,[0.1,0.5,0.95]);
    Pic{t}=pic;
end

end

