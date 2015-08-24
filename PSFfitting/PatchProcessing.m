clc;
clear;
c=0;
while c==0
    clear;
    c=0;
    m=11;
    n=11;
    mshow=10;
    nshow=16;
    [fname,pname]=uigetfile('*.tif;*.mat');
    if strcmp(fname(end-2:end),'tif')
        ImageArea=readtif([pname,fname]);
    else
        ImageArea=importdata([pname,fname]);
    end
    [Sigma,XYITlist,PeakList,Goldset]=SigmaFittingProcess(ImageArea,0,0);
    save([fname(1:end-4),'Area.mat'],'ImageArea');
    save([fname(1:end-4),'Peak.mat'],'PeakList');
    save([fname(1:end-4),'XYIT.mat'],'XYITlist');
    save([fname(1:end-4),'Gold.mat'],'Goldset');
    L=length(Goldset);
    Perm=randperm(L);
    k=1;
    ImWhole=zeros(m*mshow,n*nshow);
    PoL=zeros(0,2);
    for i=1:L
        Img=Goldset{Perm(i)};
        Po=XYITlist(Perm(i),1:2);
        if size(Img,1)==m && size(Img,2)==n
            idxx=mod(k-1,mshow);
            idxy=floor((k-1)/mshow);
            PoL=[PoL;[idxx*m+Po(1),idxy*n+Po(2)]];
            ImWhole(idxx*m+1:idxx*m+m,idxy*n+1:idxy*n+n)=Img;
            if k>mshow*nshow-1
                break;
            end
            k=k+1;
        end
    end
    axes;
    imshow(1.5*ImWhole/max(max(ImWhole)));
    hold on;
    for i=1:size(PoL,1)
        scatter(PoL(i,2),PoL(i,1),'o','b');
        hold on;
    end
    hold off;
end