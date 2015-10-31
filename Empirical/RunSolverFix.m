function [ pic,no,mv ] = RunSolverFix( bm, OptPara, StopPara, mv0 )
%   Modified in 2015.08.21, by ZHANG Haowen
%   Usage: [pic,no,mv,Time] = RunSolverFix(bm, Lambda, [sigma,bsize,psfdecay] or [sigma,bsize,bdecay,NoiseGridSize1,NoiseGridSize2], [MaxIte,Izero,Pzero], save/kframes, mv0);  
%   if OptPara=[sigma,bsize,psfdecay], then the background noise is considered to be even
%   Otherwise noise is modeled by first order b-spline base.

%%% define image size and optimization parameters
n=size(mv0.pic,1);
sigma=OptPara(1);
bsize=OptPara(2);
bdecay=OptPara(3);
s1=size(bm,1);
s2=size(bm,2);

[pic,no,mv]=EMFixSmoothMex(bm,n,sigma,StopPara,bsize,bdecay,0,[int64(1)],mv0);

end

