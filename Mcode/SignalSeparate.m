function [ Signal, Bg ] = SignalSeparate( Img, pic, no, psfno, OptPara )
%SIGNALSEPARATE Summary of this function goes here
%   Detailed explanation goes here

% Define basic parameters
s1=size(Img,1);
s2=size(Img,2);
Sigma=OptPara(1);
bsize=OptPara(2);
psfdecay=OptPara(3);
sbs=(s1+2*bsize)*(s2+2*bsize);
% Define image positions list
cx1=[0:s1-1+2*bsize]';
cx1=repmat(cx1,1,s2+2*bsize);
cx2=[0:s2-1+2*bsize];
cx2=repmat(cx2,s1+2*bsize,1);
cx=[reshape(cx1,sbs,1),reshape(cx2,sbs,1)];
% Do separation (E-step)


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
        A=sparse(psfno{i}(:,1)+1,psfno{i}(:,2)+1,psfno{i}(:,3),s1,s2);
        Bn=Bn+no(i)*full(A);
    end
    Bn=reshape(Bn,size(Bn,1)*size(Bn,2),1);
end


