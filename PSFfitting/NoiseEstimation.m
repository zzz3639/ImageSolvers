function [ NoEst, NoEstList ] = NoiseEstimation( ImageArea, PeakList )
%Modified in 2015.04, by ZHANG Haowen
%NOISEESTIMATION Summary of this function goes here
%   Estimate noise variance=NoEst^2, 
%   NoEstList is the variance list of each areas.
w =5;
wt=20; %Area of 2*w+1 in spacial and 2*wt+1 in temporal is used to do noise estimation   
M=size(ImageArea,1);
N=size(ImageArea,2);
T=size(ImageArea,3);
L=size(PeakList,1);

NoiseArea=ones(M,N,T);
NoW=zeros(L,1);
NoMean=zeros(L,1);

%Label signal areas
for l=1:L
    x=PeakList(l,1);
    y=PeakList(l,2);
    t=PeakList(l,3);
    xu=max(1,x-w);
    xd=min(M,x+w);
    yl=max(1,y-w);
    yr=min(N,y+w);
    NoiseArea(xu:xd,yl:yr,t)=0;
end

%Compute Noise Mean of each area.
for l=1:L
    x=PeakList(l,1);
    y=PeakList(l,2);
    t=PeakList(l,3);
    xu=max(1,x-w);
    xd=min(M,x+w);
    yl=max(1,y-w);
    yr=min(N,y+w);
    tb=max(1,t-wt);
    tf=min(T,t+wt);
    NoW(l)=sum(sum(sum(NoiseArea(xu:xd,yl:yr,tb:tf))));
    if NoW(l)>0
        NoMean(l)=sum(sum(sum(ImageArea(xu:xd,yl:yr,tb:tf).*NoiseArea(xu:xd,yl:yr,tb:tf))))/NoW(l);
    end
end

NoVar=zeros(L,1);
NoEstList=zeros(L,1);
for l=1:L
    x=PeakList(l,1);
    y=PeakList(l,2);
    t=PeakList(l,3);
    xu=max(1,x-w);
    xd=min(M,x+w);
    yl=max(1,y-w);
    yr=min(N,y+w);
    tb=max(1,t-wt);
    tf=min(T,t+wt);
    NoVar(l)=sum(sum(sum(NoiseArea(xu:xd,yl:yr,tb:tf).*(ImageArea(xu:xd,yl:yr,tb:tf)-NoMean(l)).*(ImageArea(xu:xd,yl:yr,tb:tf)-NoMean(l)))));
    if NoW(l)>0
        NoEstList(l)=sqrt(NoVar(l)/NoW(l));
    end
end
NoEst=sqrt(sum(NoVar)/sum(NoW));

end







