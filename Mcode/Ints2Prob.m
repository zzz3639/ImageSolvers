function [ Prob ] = Ints2Prob( Ints, Fold )
%Modified by ZHANG Haowen in 2015.09.01
%  This function converts molecule intensities to probabilities
%  Usage: [ Prob ] = Ints2Prob( Ints, Fold );     
%    Ints: a vector of Intensities;
%    Fold: a real number;
%      Here a very simple strategy is used. For Molecules brighter than
%  SoftCutoff, the probabilities are set to 1, otherwise probability equals
%  to Intensity/SoftCutoff.
%      SoftCutoff is defined by Intensity distribution: P(I>SoftCutoff)=Fold
    LogInt=log(Ints);
    M=mean(LogInt);
    S=std(LogInt);
    T=norminv(Fold);
    SoftCutoff=exp(M-S*T);
    
    Prob=Ints-(Ints>SoftCutoff).*(Ints-SoftCutoff);
    Prob=Prob/SoftCutoff;
    
end

