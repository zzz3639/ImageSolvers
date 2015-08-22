function [ pic,no,mv ] = RunSolverFix( bm, n, sigma, StopPara, optPARA, lambda, SavePara, mv0 )
%    Modified in 2015.08.21, by ZHANG Haowen
%    [pic,no,mv,Time] = RunSolverFix(bm, n, sigma, [maxite,Izero,Pzero], [bsize, bdecay], lambda, save/kframes, mv0);  


%%% define image size and optimization parameters
bsize=optPARA(1);
bdecay=optPARA(2);
s1=size(bm,1);
s2=size(bm,2);

[pic,no,mv]=EMFixSmoothMex(bm,n,sigma,StopPara,bsize,bdecay,lambda,[int64(SavePara(1))],mv0);

end

