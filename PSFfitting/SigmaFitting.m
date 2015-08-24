function [ Sigma, XYITlistFit2, PeakList3, PatchesR2 ] = SigmaFitting( ImageAreas, PeakList, alpha, NoEst0 )
%Modified in 2015.05.07, by ZHANG Haowen
%SIGMAFITTING Summary of this function goes here
%   Estimate Gaussian Sigma based on the sparse images and peak list.
%   Gaussian fitting is used here.
%   NoEst0<0 for image transfered to photon count. Otherwise NoEst0^2 should
%   be noise variance estimated outside of this program.
Maxite=1000;
M=size(ImageAreas,1);
N=size(ImageAreas,2);
T=size(ImageAreas,3);
w=5;
SignalMinFold=30; %Molecule should be brighter than 25–35 times the noise standard deviation
beta=1.96*0.6; %p=0.025, for this is a fitting based algorithm, squared error could be smaller than NoEst0.
alpha=alpha*beta;

%Do sigma free gaussian fitting and filter bad peaks
[XYITlist, Sigma] = GaussFitting( ImageAreas, PeakList, 1, w, 'free' );
L=size(XYITlist,1);
Patches=cell(L,1);
Remove=zeros(L,1);
for l=1:L
    x=PeakList(l,1);
    y=PeakList(l,2);
    t=PeakList(l,3);
    xu=max(1,x-w);
    xd=min(M,x+w);
    yl=max(1,y-w);
    yr=min(N,y+w);
    Patches{l}=ImageAreas(xu:xd,yl:yr,t);
    if NoEst0<0
        NoEst=alpha*sqrt(sum(sum(Patches{l})));
    else
        NoEst=beta*sqrt(size(Patches{l},1)*size(Patches{l},2))*NoEst0;
    end
    if SignalMinFold*NoEst/beta/sqrt(size(Patches{l},1)*size(Patches{l},2))>XYITlist(l,3) || sqrt(XYITlist(l,5))>NoEst
        Remove(l)=1;
    end
end

%Compute Sigma on remaining patches with x,y fixed
R=sum(1-Remove);
PatchesR=cell(R,1);
XYITlist2=zeros(R,size(XYITlist,2));
PeakList2=zeros(R,3);
k=1;
for l=1:L
    if Remove(l)==0
        PatchesR{k}=Patches{l};
        XYITlist2(k,:)=XYITlist(l,:);
        PeakList2(k,:)=PeakList(l,:);
        k=k+1;
    end
end
Sigma0=sqrt(sum(Sigma.*Sigma.*(1-Remove))/sum(1-Remove));
[XYITlistFit, Sigma1] = Modelfitting(PatchesR,XYITlist2,Sigma0,Maxite,'fixed');

%Filter again
LL=size(XYITlistFit,1);
Remove2=zeros(LL,1);
for l=1:LL
    if NoEst0<0
        NoEst=alpha*sqrt(sum(sum(PatchesR{l})));
    else
        NoEst=beta*sqrt(size(PatchesR{l},1)*size(PatchesR{l},2))*NoEst0;
    end
    if SignalMinFold*NoEst/beta/sqrt(size(PatchesR{l},1)*size(PatchesR{l},2))>XYITlistFit(l,3) || sqrt(XYITlistFit(l,5))>NoEst
        Remove2(l)=1;
    end
end

%Do fixed Sigma Gaussian fitting and filter again
R=sum(1-Remove2);
PatchesR2=cell(R,1);
XYITlist3=zeros(R,size(XYITlist,2));
PeakList3=zeros(R,3);
k=1;
for l=1:LL
    if Remove2(l)==0
        PatchesR2{k}=PatchesR{l};
        XYITlist3(k,:)=XYITlistFit(l,:);
        PeakList3(k,:)=PeakList2(l,:);
        k=k+1;
    end
end
[XYITlistFit2,Sigma] = Modelfitting(PatchesR2, XYITlist3, Sigma1, Maxite, 'free');

fprintf('\nPSF estimated from %d samples\n',size(XYITlistFit2,1));
fprintf('Sigma=%f\n',Sigma);

end




function [XYITlistFit, Sigma] = Modelfitting(Patches, XYITlist, Sigma0, Maxite, Opt)
%Model fitting is also solved as an MLE optimization problem. 
%EM algorithm is taken here.

%Initialize
L=length(Patches);
Boundary=4;
XYzero=1e-5;
Sigzero=1e-5;
TotalInt=0;
for i=1:L
    TotalInt=TotalInt+sum(sum(Patches{i}));
end
Izero=1e-5*TotalInt;
Nozero=1e-5*TotalInt;

XYITlistFit=XYITlist;
XYITlistFit(:,1:2)=XYITlistFit(:,1:2)+Boundary;
Sigma=Sigma0;
PatchB  =cell(size(Patches));
PatchSig=cell(size(Patches));
PatchNo =cell(size(Patches));
for i=1:L
    Img =Patches{i};
    ImgB=zeros(size(Img)+2*Boundary);
    ImgB(Boundary+1:Boundary+size(Img,1),Boundary+1:Boundary+size(Img,2))=Img;
    PatchB{i}=ImgB;
    PatchSig{i}=zeros(size(Patches{i})+2*Boundary);
    PatchNo{i} =zeros(size(Patches{i})+2*Boundary);
end

%EM Iterations
for t=1:Maxite
    %E-step
    for i=1:L
        ImgO=Patches{i};
        ImgB=PatchB{i};
        M=size(ImgO,1);
        N=size(ImgO,2);
        Mb=size(ImgB,1);
        Nb=size(ImgB,2);
        x=XYITlistFit(i,1);
        y=XYITlistFit(i,2);
        No=XYITlistFit(i,4);
        I =XYITlistFit(i,3);
        Imgx=repmat([1:Mb]',1,Nb);
        Imgy=repmat([1:Nb] ,Mb,1);
        Emol=zeros(Mb,Nb);
        Emol=I/(2*pi*Sigma*Sigma)*exp(-((Imgx-x).^2+(Imgy-y).^2)/(2*Sigma*Sigma));
        Enoi=zeros(Mb,Nb);
        Enoi=No/Mb/Nb*ones(size(Enoi));
        Esum=Enoi+Emol;
        Enoi(Boundary+1:Boundary+M,Boundary+1:Boundary+N)=ImgO.*Enoi(Boundary+1:Boundary+M,Boundary+1:Boundary+N)./Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
        Emol(Boundary+1:Boundary+M,Boundary+1:Boundary+N)=ImgO.*Emol(Boundary+1:Boundary+M,Boundary+1:Boundary+N)./Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
        PatchSig{i}=Emol;
        PatchNo{i}=Enoi;
    end
    %M-step
    XYITlistFit0=XYITlistFit;
    Sigma0=Sigma;
    SigSum=0;
    NorSum=0;
    for i=1:L
        ImgO=Patches{i};
        ImgB=PatchB{i};
        M=size(ImgO,1);
        N=size(ImgO,2);
        Mb=size(ImgB,1);
        Nb=size(ImgB,2);
        Imgx=repmat([1:Mb]',1,Nb);
        Imgy=repmat([1:Nb] ,Mb,1);
        if strcmp(Opt,'free')
            x=sum(sum(Imgx.*PatchSig{i}))/sum(sum(PatchSig{i}));
            y=sum(sum(Imgy.*PatchSig{i}))/sum(sum(PatchSig{i}));
            XYITlistFit(i,1)=x;
            XYITlistFit(i,2)=y;
        else
            x=XYITlistFit(i,1);
            y=XYITlistFit(i,2);
        end
        I=sum(sum(PatchSig{i}));
        No=sum(sum(PatchNo{i}));
        XYITlistFit(i,3)=I;
        XYITlistFit(i,4)=No;
        SigSum=SigSum+sum(sum(PatchSig{i}.*((Imgx-x).*(Imgx-x)+(Imgy-y).*(Imgy-y))));
        NorSum=NorSum+sum(sum(PatchSig{i}));
    end
    Sigma=sqrt(SigSum/NorSum/2);
    %Judge stop criterion
    if sum(abs(XYITlistFit0(:,1)-XYITlistFit(:,1)))/L<XYzero & sum(abs(XYITlistFit0(:,2)-XYITlistFit(:,2)))/L<XYzero & abs(Sigma-Sigma0)<Sigzero & ...
       sum(abs(XYITlistFit0(:,3)-XYITlistFit(:,3)))/L<Izero & sum(abs(XYITlistFit0(:,4)-XYITlistFit(:,4)))/L<Nozero
        break;
    end
end

%Compute squared error
for i=1:L
    ImgO=Patches{i};
    ImgB=PatchB{i};
    M=size(ImgO,1);
    N=size(ImgO,2);
    Mb=size(ImgB,1);
    Nb=size(ImgB,2);
    x=XYITlistFit(i,1);
    y=XYITlistFit(i,2);
    No=XYITlistFit(i,4);
    I =XYITlistFit(i,3);
    Imgx=repmat([1:Mb]',1,Nb);
    Imgy=repmat([1:Nb] ,Mb,1);
    Emol=zeros(Mb,Nb);
    Emol=I/(2*pi*Sigma*Sigma)*exp(-((Imgx-x).^2+(Imgy-y).^2)/(2*Sigma*Sigma));
    Enoi=zeros(Mb,Nb);
    Enoi=No/Mb/Nb*ones(size(Enoi));
    Esum=Enoi+Emol;
    Esum=Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
    Esum=Esum/sum(sum(Esum));
    DistDif=sum(sum((ImgO-Esum*sum(sum(ImgO))).^2));
    XYITlistFit(i,5)=DistDif;
end

XYITlistFit(:,1:2)=XYITlistFit(:,1:2)-Boundary;

end






