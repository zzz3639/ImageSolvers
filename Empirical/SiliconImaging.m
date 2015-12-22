function [ Img ] = SiliconImaging( m, s, pic, no, sigma, obnoise )
%Modified in 2015.08.21 by ZHANG Haowen
% On silicon imaging of a molecule set.
% Usage: [ Img ] = SiliconImaging( m, s, pic, no, sigma, obnoise )
%    m: number of repeats
%    s=[s1,s2], field of vision, start from 0
%    pic=[X,Y,Intensity], molecule list
%    no: noise, count by photon number
%    sigma: PSF width
%    obnoise: Additional noise, for each pixel, signal=photoncount+sqrt(photoncount)*rand()*obnoise

%%% initialize
rng('shuffle');
n=size(pic,1);
s1=s(1);
s2=s(2);
% define image positions list
cx1=[0:s1-1]';
cx1=repmat(cx1,1,s2);
cx2=[0:s2-1];
cx2=repmat(cx2,s1,1);
cx=[reshape(cx2,s1*s2,1),reshape(cx1,s1*s2,1)];
% define output variable
Img=cell(m,1);

%%% generate poisson images
for i=1:m
    img=zeros(s1*s2,1);
    for j=1:n
        img=img+[mnrnd(round(pic(j,3)),PSF(pic(j,1:2),cx,sigma))]';
    end
    img=img+[Noise(s1,s2,round(no))]';
    Img{i}=reshape(img,s1,s2);
end

%%% add additional noise
for i=1:m
    img=Img{i};
    imgb=zeros(size(img));
    imgb=img+obnoise*randn(s1,s2).*sqrt(img);
    imgb=imgb-(imgb<0).*imgb;
    Img{i}=imgb;
end

end

function b2=Noise(s1,s2,PhotonNum)
    b2=mnrnd(PhotonNum,repmat(1/s1/s2,s1*s2,1));
end

%PSF function
function V=PSF(mu,cx,sig)
    V=normpdf(cx(:,1),mu(1),sig).*normpdf(cx(:,2),mu(2),sig);
    V=V/sum(V);
end


