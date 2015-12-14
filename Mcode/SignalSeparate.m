function [ Signal, Bg, StandardPSF, ExpectedSignal, ExpectedBg ] = SignalSeparate( Img, pic, no, psfno, OptPara )
%SIGNALSEPARATE Summary of this function goes here
%   Detailed explanation goes here

% Define basic parameters
s1=size(Img,1);
s2=size(Img,2);
Sigma=OptPara(1);
bsize=OptPara(2);
psfdecay=OptPara(3);
sb1=(s1+2*bsize);
sb2=(s2+2*bsize);
sbs=(s1+2*bsize)*(s2+2*bsize);
n=size(pic,1);
B=reshape(Img,s1*s2,1);

pic=[pic(:,1),pic(:,2),pic(:,3)];

cx1=[0:s1-1+2*bsize]';
cx1=repmat(cx1,1,s2+2*bsize);
cx2=[0:s2-1+2*bsize];
cx2=repmat(cx2,s1+2*bsize,1);
cx=[reshape(cx1,sbs,1),reshape(cx2,sbs,1)];

if size(psfno,1)==0
    tempno=zeros(sbs,3);
    tempno(:,1:2)=cx;
    tempno(:,3)=1/sbs;
    psfno=cell(1,1);
    psfno{1}=tempno;
end
% Do separation (E-step)
A =PSF(0,s1-1,0,s2-1,pic(:,1:2),Sigma,psfdecay);
Bn=PSFno(psfno,no,sb1,sb2);
BnFull=reshape(Bn,sb1,sb2);
Bn=reshape(BnFull(bsize:bsize+s1-1,bsize:bsize+s2-1),s1*s2,1);
SA=A*pic(:,3)+Bn;
W=spdiags(B./SA(:,1),0,s1*s2,s1*s2) * A * spdiags(pic(:,3),0,n,n);
Wn=Bn./SA;
Signal=W;
Bg=Wn.*B;
ExpectedSignal = A*spdiags(pic(:,3),0,n,n);
ExpectedBg = Bn;
StandardPSF=PSF(0,s1-1,0,s2-1,pic(:,1:2),Sigma,psfdecay);
end

function Mu=EMPositionStep(Img,pic,no,psfno,OptPara,cx)
%Define parameters
    cx=[cx(:,2),cx(:,1)];
    s1=size(Img,1);
    s2=size(Img,2);
    Sigma=OptPara(1);
    bsize=OptPara(2);
    psfdecay=OptPara(3);
    sb1=(s1+2*bsize);
    sb2=(s2+2*bsize);
    sbs=(s1+2*bsize)*(s2+2*bsize);
    n=size(pic,1);
    Zero=1e-12;
%Put the image into a frame with margin
    bm=zeros(sbs,1);
    P1=zeros(s1+2*bsize,s2+2*bsize);
    P1(bsize+1:bsize+s1,bsize+1:bsize+s2)=Img;
    bm=reshape(P1,(s1+2*bsize)*(s2+2*bsize),1);
%E-step
    A =PSF(0,s1-1+2*bsize,0,s2-1+2*bsize,pic(:,1:2)+bsize,Sigma,psfdecay);
    Bn=PSFno(psfno,no,sb1,sb2);
    SA=A*pic(:,3)+Bn;
    W= spdiags(ones(sbs,1)./SA,0,sbs,sbs) * A * spdiags(pic(:,3),0,n,n);
    Wn=Bn./SA;
    Judge=(cx(:,1)>=bsize).*(cx(:,1)<bsize+s2).*(cx(:,2)>=bsize).*(cx(:,2)<bsize+s1);
    bbm=SA.*(1-Judge)+bm.*Judge;        
%M-step
    Mu=zeros(n,2);
    %Update first axes
    T=(W)'*sparse(bbm.*cx(:,1));
    O=(W)'*sparse(bbm);
    for j=1:n
        if O(j,1)>Zero
            P(j)=O(j,1)^(-1)*T(j,1);
        end
    end
    Mu(:,1)=P'-bsize;
    %Update second axes
    T=(W)'*sparse(bbm.*cx(:,2));
    O=(W)'*sparse(bbm);
    for j=1:n
        if O(j,1)>Zero
            P(j)=O(j,1)^(-1)*T(j,:);
        end
    end
    Mu(:,2)=P'-bsize;
end

function Bp=PSF(L1, U1, L2, U2, Mu, Sigma,bdecay)
    n1=U1-L1+1;
    n2=U2-L2+1;
    nspots=size(Mu,1);
    %Bp=sparse(n1*n2,1);
    I=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    J=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    V=zeros(nspots*4*(bdecay+1)*(bdecay+1),1);
    num=0;
    for k=1:nspots
        m1=floor(Mu(k,1));
        m2=floor(Mu(k,2));
        l1=max(m1-bdecay,L1);
        u1=min(m1+bdecay+1,U1);
        l2=max(m2-bdecay,L2);
        u2=min(m2+bdecay+1,U2);
        t1=u1-l1+1;
        t2=u2-l2+1;
        if t1<0
            t1=0;
        end
        if t2<0
            t2=0;
        end
    %for i=l1:u1
    %    Bp((i-L1)*n2+l2-L2+1:(i-L1)*n2+u2-L2+1,1)=1/(2*pi*sigma*sigma)*exp((-(i-mu(1))^2-([l2:u2]'-mu(2)).^2)/(2*sigma*sigma));
    %end
        Bp1=zeros(1,u1-l1+1);
        Bp1(1,1:u1-l1+1)=exp(-(([l1:u1]-Mu(k,1)).^2)/(2*Sigma*Sigma));
        Bp2=zeros(u2-l2+1,1);
        Bp2(1:u2-l2+1,1)=exp(-(([l2:u2]'-Mu(k,2)).^2)/(2*Sigma*Sigma));
        V(num+1:num+t1*t2)=reshape([(Bp2*Bp1)/(2*pi*Sigma*Sigma)]',t1*t2,1);
        II=(repmat([l1:u1]',1,t2)-L1)*n2+repmat([l2:u2],t1,1)-L2+1;
        I(num+1:num+t1*t2)=reshape(II,t1*t2,1);
        J(num+1:num+t1*t2)=k;
        num=num+t1*t2;
    %ones(n1,1)*ones(1,n2);
    end
    Bp=sparse(I(1:num),J(1:num),V(1:num),n1*n2,nspots);
end

function Bn=PSFno(psfno,no,sb1,sb2)
    Bn=zeros(sb1,sb2);
    for i=1:length(psfno)
        A=sparse(psfno{i}(:,1)+1,psfno{i}(:,2)+1,psfno{i}(:,3),sb1,sb2);
        Bn=Bn+no(i)*full(A);
    end
    Bn=reshape(Bn,size(Bn,1)*size(Bn,2),1);
end


