psfdecay=5;
s1=17;
s2=17;
Mu=Pic(j,1:2);
U1=floor(Mu(2))+1-psfdecay; U1=max(1,U1);
L1=floor(Mu(2))+2+psfdecay; L1=min(s1,L1);
U2=floor(Mu(1))+1-psfdecay; U2=max(1,U2);
L2=floor(Mu(1))+2+psfdecay; L2=min(s2,L2);
T1=reshape(P1,s1,s2);
T2=reshape(P2,s1,s2);
PlotShow=zeros(L1-U1+1,2*(L2-U2+1));
PlotShow=[T1(U1:L1,U2:L2),T2(U1:L1,U2:L2)];
imshow(full(PlotShow));