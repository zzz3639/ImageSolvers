function Res=removezeros(Pic,fold)
% Modified in 2015.06.20 by ZHANG Haowen
% Usage: Res=removezeros(Pic,fold)
%   Pic: Molecule list, [x,y,intensity]
%   fold: Molecule intensities are assumed to follow log normal distribution. 
%         1-fold is the p-value in one side statistic test
    Int=Pic(:,3);
    LogInt=log(Int);
    M=mean(LogInt);
    S=std(LogInt);
    T=norminv(fold);
    Cutoff=M-S*T;
    Remain=(LogInt>Cutoff);
    Res=Pic(Remain,:);
end