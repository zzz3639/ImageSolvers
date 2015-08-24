function [ FPeakList ] = PeakFilter( ImageArea, PeakList )
%PEAKFILTER Summary of this function goes here
%   This function remove peaks likely to be false.
sz=3;
w=5;
k=1;
M=size(ImageArea,1);
N=size(ImageArea,2);
T=size(ImageArea,3);

% Average of the adjacent area (2*w+1)*(2*w+1) area is a upper estimation of the local noise.
% The second largest pixel around a peak is required to be >noise+th. The several pixels around the true peak should also be large.
% For poisson noise, variance=mean, th=1.96*sqrt(no).
avmask=zeros(M,N,T);
PicThis=zeros(M,N);
avpic=ones(M,N);
avker=ones(2*w+1,2*w+1);
% Compute centered average
for t=1:T
    PicThis=ImageArea(:,:,t);
    PicDiv=conv2(avpic,avker);
    PicSum=conv2(PicThis,avker);
    PicAve=PicSum./PicDiv;
    PicAve=PicAve(w+1:end-w,w+1:end-w);
    avmask(:,:,t)=PicAve;
end

%Compare peaks with local averages
FPeakList=zeros(0,3);
L=size(PeakList,1);
for i=1:L
    x=PeakList(i,1);
    y=PeakList(i,2);
    t=PeakList(i,3);
    S=reshape(ImageArea(x-k:x+k,y-k:y+k,t),(2*k+1)*(2*k+1),1);
    S=sort(S);
    if S(end-1)>avmask(x,y,t)+1.96*sqrt(avmask(x,y,t));
        FPeakList=[FPeakList;[x,y,t]];
    else
    end
end

% Nearby peaks are considered as one
% Locate the Peaks in this image area 
PeakPic=zeros(M,N,T);
L=size(FPeakList,1);
for i=1:L
    x=FPeakList(i,1);
    y=FPeakList(i,2);
    t=FPeakList(i,3);
    PeakPic(x,y,t)=ImageArea(x,y,t);
end

% Remove the peaks with brighter peaks nearby
for i=1:L
    x =FPeakList(i,1);
    y =FPeakList(i,2);
    t =FPeakList(i,3);
    xu=max(x-sz,1);
    xd=min(x+sz,M);
    yl=max(y-sz,1);
    yr=min(y+sz,N);
    Peakthis=PeakPic(xu:xd,yl:yr,t);
    Peakthis=reshape(Peakthis,size(Peakthis,1)*size(Peakthis,2),1);
    Peakthis=sort(Peakthis);
    if PeakPic(x,y,t)<Peakthis(end) || Peakthis(end)==Peakthis(end-1)
        PeakPic(x,y,t)=0;
    end
end

FPeakList=zeros(0,3);
for t=1:T
    [I,J]=find(PeakPic(:,:,t));
    FPeakList=[FPeakList;[I,J,t*ones(length(I),1)]];
end

end







