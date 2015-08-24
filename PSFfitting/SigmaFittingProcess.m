function [ Sigma, XYITlist, PeakList, Goldset ] = SigmaFittingProcess( ImageArea, PorN, Seed )
%Modified in 2015.04, by ZHANG Haowen
%SIGMAF Summary of this function goes here
%   Fitting Gaussian Sigma from image series.
%   If PorN==1, that means this image signal is transfered to photon count,
% then background noise is assumed to be 1.5*Possion.
%   If PorN==0, variance of background noise is estimated from the Image series
%   Seed is random seed
rng(Seed);
PeakList=PeakCalling(ImageArea, 0);
PeakList=PeakFilter(ImageArea, PeakList);

if PorN==1
    NoEst=sqrt(1.5*sum(sum(sum(ImageArea)))/size(ImageArea,1)/size(ImageArea,2)/size(ImageArea,3));
else
    NoEst=NoiseEstimation(ImageArea,PeakList);
end

[Sigma,XYITlist,PeakList,Goldset] = SigmaFitting(ImageArea,PeakList,1,NoEst);

end

