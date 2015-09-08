function [ pic,no,mv ] = RunSolverFix( bm, OptPara, StopPara, mv0 )
%   Modified in 2015.08.21, by ZHANG Haowen
%   Usage: [pic,no,mv,Time] = RunSolverFix(bm, Lambda, [sigma, bsize, bdecay], [MaxIte,Izero,Pzero], save/kframes, mv0);  


%%% define image size and optimization parameters
n=size(mv0.pic,1);
sigma=OptPara(1);
bsize=OptPara(2);
bdecay=OptPara(3);
s1=size(bm,1);
s2=size(bm,2);

[pic,no,mv]=EMFixSmoothMex(bm,n,sigma,StopPara,bsize,bdecay,0,[int64(1)],mv0);

end

