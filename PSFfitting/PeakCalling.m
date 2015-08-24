function [ PeakList, Masks ] = PeakCalling( ImageArea , th)
%Modified in 2015.04.28, by ZHANG Haowen
%PEAKCALLING Summary of this function goes here
%   This function aims to find peaks in image series.
%   Peaks found here undergoes a preliminary filter. 
%   Pixels>(max-average+)/2+average are considered as signal
%   then a svm classifier is trained to obtain all the signal pixels.

Zero=1e-6;
%medium filter radius
%notice that peak are required to locate k+1 pixels away from the boundary.
k=1;

M=size(ImageArea,1);
N=size(ImageArea,2);
T=size(ImageArea,3);

Masks=zeros(M,N,T);
Mask  =zeros(M,N);
Peaks =zeros(M,N);
PeaksU=zeros(M,N);
PeaksD=zeros(M,N);
PeaksL=zeros(M,N);
PeaksR=zeros(M,N);
PeaksUL=zeros(M,N);
PeaksUR=zeros(M,N);
PeaksDL=zeros(M,N);
PeaksDR=zeros(M,N);
PeaksB=zeros(M,N);

PeaksImg=zeros(M,N,T);
LabelImg=zeros(M,N,T);
SVMImg=zeros(M,N,T);

TrainSery=zeros((M-2)*(N-2),(2*k+1)*(2*k+1),T);

%Find peaks and find positive labels.
for t=1:T
    %smooth by median filter
    [Temppic,PixelData]=Median(ImageArea(:,:,t),th,2*k+1);
    TrainSery(:,:,t)=PixelData;
    Mask(k+1:M-k,k+1:N-k)   =ImageArea(k+1:M-k,k+1:N-k,t);
    Masks(k+1:M-k,k+1:N-k,t)=ImageArea(k+1:M-k,k+1:N-k,t);
    %find the peaks
    for i=k+2:M-k-1
        PeaksL(i,k+2:N-k-1)=(Mask(i,k+2:N-k-1)>=Mask(i,k+1:N-k-2));
        PeaksR(i,k+2:N-k-1)=(Mask(i,k+2:N-k-1)>=Mask(i,k+3:N-k));
    end
    for j=k+2:N-k-1
        PeaksU(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+1:M-k-2,j));
        PeaksUL(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+1:M-k-2,j-1));
        PeaksUR(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+1:M-k-2,j+1));
        PeaksD(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+3:M-k,j));
        PeaksDL(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+3:M-k,j-1));
        PeaksDR(k+2:M-k-1,j)=(Mask(k+2:M-k-1,j)>=Mask(k+3:M-k,j+1));
    end
    PeaksB=Mask>Zero;
    Temppic=PeaksB.*PeaksU.*PeaksD.*PeaksL.*PeaksR.*PeaksUL.*PeaksUR.*PeaksDL.*PeaksDR;
    PeaksImg(:,:,t)=Temppic;
    %find positive labels
    Ave=mean(mean(ImageArea(:,:,t)));
    Max=max(max(max(ImageArea(:,:,:))));
    ThImg=Ave+(Max-Ave)/2;
    LabelImg(:,:,t)=(ImageArea(:,:,t)>ThImg).*Temppic;
end

%Train SVM classifier as a preliminary filter. It can be proved that about
%1/k/k proportion of pixels are local maximum, a preliminary filter significantly improve the speed.
pixelsum=(M-2*k)*(N-2*k)*T;
Plabelsum=sum(sum(sum(LabelImg)));
Xposi=zeros(Plabelsum,(2*k+1)*(2*k+1));
t=ceil(pixelsum/(pixelsum-Plabelsum)*Plabelsum);
NegaRand=[randi(M-2*k,t,1)+k,randi(N-2*k,t,1)+k,randi(T,t,1)];
Nlabelsum=0;
for t=1:T
    [I]=find(NegaRand(:,3)==t);
    X=NegaRand(I,1);
    Y=NegaRand(I,2);
    Img=LabelImg(:,:,t);
    Nlabelsum=Nlabelsum+sum(1-Img((Y-1)*M+X));
end
Xnega=zeros(Nlabelsum,(2*k+1)*(2*k+1));
lp=0;
ln=0;
for t=1:T
    Img=LabelImg(:,:,t);
    [Ip]=find(Img);
    Indn=find(NegaRand(:,3)==t);
    X=NegaRand(Indn,1);
    Y=NegaRand(Indn,2);
    [Indn2]=find(1-Img((Y-1)*M+X));
    In=X(Indn2,:);
    Jn=Y(Indn2,:);
    In=(Jn-1)*M+In;
    Np=length(Ip);
    Nn=length(In);
    Img=ImageArea(:,:,t);
    for i=[-k:k]
        Xposi(lp+1:lp+Np,(i+k)*(2*k+1)+1:(i+k+1)*(2*k+1))=(Img(repmat((Ip+i*M)',2*k+1,1)+repmat([-k:k]',1,length(Ip))))';
        Xnega(ln+1:ln+Nn,(i+k)*(2*k+1)+1:(i+k+1)*(2*k+1))=(Img(repmat((In+i*M)',2*k+1,1)+repmat([-k:k]',1,length(In))))';
    end
    lp=lp+Np;
    ln=ln+Nn;
end

%training SVM using libsvm package
%ModelSVM=svmtrain([zeros(Nlabelsum,1);ones(Plabelsum,1)],[Xnega;Xposi],'-t 0 -c 1');
%svmW=((ModelSVM.sv_coef)'*ModelSVM.SVs)';
%svmB=ModelSVM.rho;

%sort the dimensions
Xposi=sort(Xposi,2);
Xnega=sort(Xnega,2);
%training L2-regularized logistic regression using liblinear.
Model=train([zeros(Nlabelsum,1);ones(Plabelsum,1)],sparse([Xnega;Xposi]),'-s 0 -c 100 -e 0.00000001 -B 1');
W=Model.w;
svmW=(W(1:end-1))';
svmB=W(end);

%training logistic regression using matlab mnrfit
%Model=mnrfit([Xnega;Xposi],[zeros(Nlabelsum,1),ones(Nlabelsum,1);ones(Plabelsum,1),zeros(Plabelsum,1)]);
%svmW=Model(2:end,:);
%svmB=Model(1,1);

if sum(Xposi*svmW-svmB)<0
    svmW=-svmW;
    svmB=-svmB;
end

%classify the image series
for t=1:T
    PixelData=sort(TrainSery(:,:,t),2);
    ThisLabel=((PixelData*svmW+svmB)>0);
    SVMImg(k+1:M-k,k+1:N-k,t)=reshape(ThisLabel,M-2*k,N-2*k);
end

PeakList=zeros(0,3);
for t=1:T
    Temppic=PeaksImg(:,:,t).*SVMImg(:,:,t);
    Mask=Masks(:,:,t);
    [I,J]=find(Temppic);

    for i=1:length(I)
        x=I(i);
        y=J(i);
        if Mask(x,y)>=Mask(x,y-1) && Mask(x,y)>=Mask(x-1,y) && Mask(x,y)>=Mask(x+1,y) && Mask(x,y)>=Mask(x,y+1)  && ...
           Mask(x,y)>=Mask(x+1,y+1) && Mask(x,y)>=Mask(x-1,y-1) && Mask(x,y)>=Mask(x+1,y-1) && Mask(x,y)>=Mask(x-1,y+1) && ...
           Mask(x,y)>Zero    %It is possible to add more conditions here
            PeakList=[PeakList;[x,y,t]];
        else
            Temppic(x,y)=0;
        end
    end
    
end

end


function [MF,TrainData]=Median(img, th, k)
m=size(img,1);
n=size(img,2);
T1=zeros((m+1-k)*k,n);
for i=1:k
    T1(i:k:end,:)=img(i:m+i-k,:);
end
T2=zeros((m+1-k)*k,(n+1-k)*k);
for i=1:k
    T2(:,i:k:end)=T1(:,i:n+i-k);
end

T1=zeros((m+1-k)*(n+1-k)*k,k);
for i=1:n+1-k
    T1((i-1)*k*(m+1-k)+1:i*k*(m+1-k),:)=T2(:,(i-1)*k+1:i*k);
end
T2=reshape(T1',k*k,(m+1-k)*(n+1-k));
TrainData=T2'; %For SVM training, nothing to do with medium filter.
mf=median(T2,1);
MF=reshape(mf,m+1-k,n+1-k);

MF=MF-th;
MF=MF-(MF<0).*MF;

end

