function [ imgpsfno ] = DrawPSFno( psfno, no, sz )
%Visualilization of noise pattern
%   Usage: [imgpsfno] = DrawPSFno( psfno, no, sz )
s1=sz(1);
s2=sz(2);

imgpsfno=zeros(s1,s2);
for i=1:length(psfno)
    A=sparse(psfno{i}(:,1)+1,psfno{i}(:,2)+1,psfno{i}(:,3),s1,s2);
    imgpsfno=imgpsfno+no(i)*full(A);
end

end

