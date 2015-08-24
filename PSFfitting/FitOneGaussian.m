function [ X,Y,Int,No,Sig,DistDif ] = FitOneGaussian( Patch, Sigma, X0, Y0, Maxite, Opt )
%Modified in 2015.04, by ZHANG Haowen
%FIXONE Summary of this function goes here
%  X,Y: Position.  Int: intensity.  No: noise.   Sig: Sigma.
%  DistDif: squared error.
%  Expected distribution: Int*normpdf(X,Y,Sig)+No/sizeof(Image Patch)
%  Gaussian MLE fitting, Fitting one spot by using EM algorithm
%  Signal is is assumed to be psf+b.
%  Opt='fixed' or 'free', 'fitted' for fitting the model with fitted sigma, 
%  'free' if we treat the sigma as free parameter and use the input as initial point.
%  MLE estimator is non-convex, thus starting point may alter the result.
%  The optimization is performed from 'Trial' number of starting points 

Trial=5;
Boundary=4;
XYzero=1e-5;
Sigzero=1e-3;
Izero=1e-5*sum(sum(Patch));
Nozero=1e-5*sum(sum(Patch));
SigmaLB=0.7; %Lower bound of sigma. Such Gaussian peak is too sharp thus considered as noise.

M =size(Patch,1);
N =size(Patch,2);
Mb=M+Boundary*2;
Nb=N+Boundary*2;
Patchb=zeros(Mb,Nb);
Imgx=repmat([1:Mb]',1,Nb);
Imgy=repmat([1:Nb] ,Mb,1);
Patchb(Boundary+1:Boundary+M,Boundary+1:Boundary+N)=Patch;
TotalPhoton=sum(sum(Patch))*Mb*Nb/M/N;
Emol=zeros(Mb,Nb);
Enoi=zeros(Mb,Nb);

Lmax=-inf;
 for k=1:Trial
    SigEnd=0;
    no=rand()*TotalPhoton;
    I=TotalPhoton-no;
    x=randn()+X0+Boundary;
    y=randn()+Y0+Boundary;
    sig=Sigma;
    for t=1:Maxite
        %E-step
        Enoi=no/Mb/Nb*ones(size(Enoi));
        Emol=I/(2*pi*sig*sig)*exp(-((Imgx-x).^2+(Imgy-y).^2)/(2*sig*sig));
        Esum=Enoi+Emol;
        Enoi(Boundary+1:Boundary+M,Boundary+1:Boundary+N)=Patch.*Enoi(Boundary+1:Boundary+M,Boundary+1:Boundary+N)./Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
        Emol(Boundary+1:Boundary+M,Boundary+1:Boundary+N)=Patch.*Emol(Boundary+1:Boundary+M,Boundary+1:Boundary+N)./Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
        %M-step
        It=I; xt=x; yt=y; noT=no; sigt=sig;
        no=sum(sum(Enoi));
        I =sum(sum(Emol));
        x=sum(sum(Emol.*Imgx))/sum(sum(Emol));
        y=sum(sum(Emol.*Imgy))/sum(sum(Emol));
        if strcmp(Opt,'fixed')
            sig;
        else
            sig=sqrt(sum(sum(Emol.*((Imgx-x).*(Imgx-x)+(Imgy-y).*(Imgy-y))))/sum(sum(Emol))/2);
        end
        %judge stop condition
        if sig<0.7
            SigEnd=1;
            break;
        end
        if abs(x-xt)<XYzero & abs(y-yt)<XYzero & abs(I-It)<Izero & abs(no-noT)<Nozero & abs(sig-sigt)<Sigzero
            break;
        end
    end
    if SigEnd~=1
        Enoi=no/Mb/Nb*ones(size(Enoi));
        Emol=I/(2*pi*sig*sig)*exp(-((Imgx-x).^2+(Imgy-y).^2)/(2*sig*sig));
        Esum=Enoi+Emol;
        Esum=Esum(Boundary+1:Boundary+M,Boundary+1:Boundary+N);
        Esum=Esum/(sum(sum(Esum)));
        L=sum(sum(Patch.*log(Esum)));
        if L>Lmax
            Lmax=L;
            X=x-Boundary;
            Y=y-Boundary;
            Int=I;
            No=no;
            Sig=sig;
            DistDif=sum(sum((Patch-Esum*sum(sum(Patch))).^2));
        end
    else
        continue;
    end
end
if ~exist('X')
    X=0;
    Y=0;
    Int=0;
    No=0;
    Sig=0;
    DistDif=inf;
end

end




