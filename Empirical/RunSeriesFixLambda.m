function [ Pic ] = RunSeriesFixLambda( Sery, sig, optPARA, Lambda )
%RUNSERIES Summary of this function goes here
%   Detailed explanation goes here
N=size(Sery,3);
Izero=1e-5;
Pzero=1e-4;
MaxIteIni=500;
MaxIte0=5000;
Izero0=1e-6;
Pzero0=1e-5;
bsize=optPARA(1);
bdecay=optPARA(2);
if length(optPARA)==2
    EvenNoise=1;
else
    EvenNoise=0;
    NoiseGridSize1=optPARA(3);
    NoiseGridSize2=optPARA(4);
end
Pic=cell(N,1);

for t=1:N
    Img=Sery(:,:,t);
    [pic,no]=FastSolver( Img, floor(size(Img,1)*size(Img,2)/7/7/sig/sig*20), Lambda, [sig,bsize,bdecay,NoiseGridSize1,NoiseGridSize2], [MaxIteIni,MaxIte0,Izero0*sum(sum(Img)),Pzero0], [0.1,0.5] );
    pic=PostRun(pic,[1.5/8,0.5,0.95]);
    
    Pic{t}=pic;
end

end

