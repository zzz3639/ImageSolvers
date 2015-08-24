function [ XYITlist, Sigma ] = GaussFitting( ImageAreas, PeakList, Sigma0, w, Opt )
%Modified in 2015.04, by ZHANG Haowen
%GAUSSFITTINGFIX Summary of this function goes here
%   Fitting the precious positions and intensities based on the peak list.
%   Opt='fixed' for sigma fitted to be Sigma0, Opt='free' if sigma is
%   treated as free variable.

MaxIte=1000;
M=size(ImageAreas,1);
N=size(ImageAreas,2);
%w=5;

L=size(PeakList,1);
XYITlist=zeros(L,5);
Sigma=zeros(L,1);
fprintf('\nGaussian Fitting Started. %d candidates in total\n',L);
for l=1:L
    if mod(l,500)==0
        l
    end
    x=PeakList(l,1);
    y=PeakList(l,2);
    t=PeakList(l,3);
    xu=max(1,x-w);
    xd=min(M,x+w);
    yl=max(1,y-w);
    yr=min(N,y+w);
    Patch=ImageAreas(xu:xd,yl:yr,t);
    [xfit, yfit, Intfit, Nofit, Sigfit, DistDif]=FitOneGaussian(Patch,Sigma0,x-xu+1,y-yl+1,MaxIte,Opt);
    XYITlist(l,:)=[xfit,yfit,Intfit,Nofit,DistDif];
    Sigma(l)=Sigfit;
end

end

