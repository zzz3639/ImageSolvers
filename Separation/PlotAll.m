psfdecay=4;
s1=17;
s2=17;
sn1=10;
sn2=10;
PlotShowAll=zeros(3*(psfdecay+1)*sn1,4*(psfdecay+1)*sn2);
kk=0;
for j=1:size(Signal,2)
    kk1=floor(kk/sn1);
    kk2=mod(kk,sn1);
    P1=Signal(:,j);
    P1=P1/sum(P1);
    P2=StandardPSF(:,j);
    P2=P2/sum(P2);
    Mu=Pic(j,1:2);
    U1=floor(Mu(2))+1-psfdecay; U1=max(1,U1);
    L1=floor(Mu(2))+2+psfdecay; L1=min(s1,L1);
    U2=floor(Mu(1))+1-psfdecay; U2=max(1,U2);
    L2=floor(Mu(1))+2+psfdecay; L2=min(s2,L2);
    T1=reshape(P1,s1,s2);
    T2=reshape(P2,s1,s2);
    T3=Img;
    PlotShowThis=zeros(2*(psfdecay+1),6*(psfdecay+1));
    PlotShowThis(1:L1-U1+1,1:L2-U2+1)=T1(U1:L1,U2:L2)*1.5;
    PlotShowThis(1:L1-U1+1,2*(psfdecay+1)+1:2*(psfdecay+1)+L2-U2+1)=T2(U1:L1,U2:L2)*1.5;
    PlotShowThis(1:L1-U1+1,4*(psfdecay+1)+1:4*(psfdecay+1)+L2-U2+1)=T3(U1:L1,U2:L2)/max(max(T3(U1:L1,U2:L2)))/1.5;
    PlotShowAll(kk1*2*(psfdecay+1)+1:(kk1+1)*2*(psfdecay+1),kk2*6*(psfdecay+1)+1:(kk2+1)*6*(psfdecay+1))=PlotShowThis;
    kk=kk+1;
end
imshow(PlotShowAll*2);