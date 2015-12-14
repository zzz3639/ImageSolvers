function [ Loss, Losses ] = PrecisionTrue( RLassoAns, TrueAns )
%PRECISIONTRUE Summary of this function goes here
%   Detailed explanation goes here
L=length(RLassoAns);
Loss=zeros(L,1);
Losses=cell(L,1);
for i=1:L
    Pic=RLassoAns{i}.pic;
    [AlignThis]=MolListAlign(Pic, TrueAns);
    Prob=Ints2Prob(AlignThis(:,3),0.75);
    E=sum(Prob.*AlignThis(:,4))/sum(Prob);
    Losses{i}=[Prob,AlignThis(:,4)];
    %E=mean(AlignThis(:,4));
    Loss(i)=E;
end

end

