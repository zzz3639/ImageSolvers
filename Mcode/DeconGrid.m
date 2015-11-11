function [ HighImg, mv ] = DeconGrid( Img, OptPara, Lambda, Grid, MaxIte, EvenNoise, MvRate, NoiseGridSize )
% Do deconvolution using SparseGMM
%   Usage: [ HighImg, mv ] = DeconGrid( Img, OptPara, Lambda, Grid, MaxIte, EvenNoise, MvRate, NoiseGridSize )
%      OptPara=[Sigma, boundarysize, psfdecay] 
s1=size(Img,1);
s2=size(Img,2);

n1=s1*Grid;
n2=s2*Grid;

if EvenNoise==0
    NoiseGridSize1=NoiseGridSize(1);
    NoiseGridSize2=NoiseGridSize(2);
end

Sigma=OptPara(1);
bsize=OptPara(2);
psfdecay=OptPara(3);

cx1temp=([0:n1-1]+0.5)/2;
cx1temp=repmat(cx1temp',1,n2);
cx2temp=([0:n2-1]+0.5)/2;
cx2temp=repmat(cx2temp, n1,1);
cx=[reshape(cx2temp,n1*n2,1),reshape(cx1temp,n1*n2,1)];

Int=0.5*ones(n1*n2,1)*mean(mean(Img));
if EvenNoise==1
    no0=0.5*sum(sum(Img));
else
    numg1=floor((s1-1-floor(NoiseGridSize1/2))/NoiseGridSize1)+1;
    numg2=floor((s2-1-floor(NoiseGridSize2/2))/NoiseGridSize2)+1;
    NumNo=numg1*numg2;
    no0=0.5/NumNo*ones(NumNo,1)*sum(sum(Img));
end

mv0.pic=[cx,Int];
mv0.no=no0;

if EvenNoise==1
    [dpic,no,mvpic]=EMFixSmoothMex(Img,n1*n2,Sigma,[MaxIte,1e-12,1e-12],bsize,psfdecay,Lambda,[int64(MvRate)],mv0);
else
    [dpic,no,mvpic]=EMSparseUnevenFixMex(Img,n1*n2,Sigma,[MaxIte,1e-12,1e-12],bsize,psfdecay,[NoiseGridSize1,NoiseGridSize2],Lambda,[int64(MvRate)],mv0);
end
mv=zeros(n1,n2,mvpic.N);
for i=1:mvpic.N
    mv(:,:,i)=reshape(mvpic.pic(:,3,i),n1,n2);
end

HighImg=reshape(dpic(:,3),n1,n2);

end

