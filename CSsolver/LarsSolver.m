function [ X, LassoPath ] = LarsSolver( A, b, c, Lambda, Type, X0 )
%Solves this Problem by Least angle regression: 
%    Type=0: 0.5*|(Ax-b)|^2+Lambda*c'*x  s.t. x>=0
%    Type=1: c'*x  s.t. |(Ax-b)|^2 <= Lambda^2, && x>=0

% A,b,c and Lambda are assumed to be non-negative.
%   
% Initial point x0 could be given outside. x0 should be the solution of the
% following optimization problem:
%   minimize: 0.5*|(Ax-b)|^2  s.t. x>=0, x_i=0 if c_i~=0
% If x0 isn't given, it is computed by function quadprog().

Zero=1e-12; % minimal value of A, below which regarded as zero
Czero=1e-15; % minimal value of c
Qzero=1e-10; % value to judge active elements in x0
m=size(A,1);
n=size(A,2);
ActiveNew=zeros(n,1);
ActiveNew=ActiveNew>0.5;
XNew=zeros(n,1);

%%% Initialize
Atb=A'*b;
Active=(c<Czero);
Idxthis=find(Active);
Athis=A(:,Active);
if length(Idxthis)>0
    if exist('X0')
    else
        H=Athis'*Athis;
        F=-Athis'*b;
        Ap=-eye(length(Idxthis));
        Bp=0;
        Xtemp=quadprog(H,F,Ap,Bp);
        X0=zeros(size(XNew));
        X0(Idxthis,:)=Xtemp;
    end
    ActiveNew(X0>=Zero,:)=1;
    XNew(X0>=Qzero,:)=X0(X0>=Qzero,:);
    AtAX=A'*(Athis*XNew(Idxthis,:));
    ValueLambda=(Atb(~Active,:)-AtAX(~Active,:))./c(~Active,:);
    u=max(ValueLambda);
    LNew=max(u,0)+1;
else
    ValueLambda=Atb./c;
    [u]=max(ValueLambda);
    ActiveNew=(ValueLambda==u);
    LNew=u;
    XNew=zeros(size(XNew));
end


L=LNew;
X=XNew;
Active=ActiveNew;
Node.X=X;
Node.L=L;
Node.Active=Active;
LassoPath=cell(1,1);
LassoPath{1}=Node;

%%% Run lars iterations
while 1
    if Type==0
        if LNew<Lambda
            Athis=A(:,Active);
            Idxthis=find(Active);
            Inv=(Athis'*Athis)^-1;
            SumInv=Inv*c(Active,:);
            Atb=Athis'*b;
            SinvAtb=Inv*Atb;
            % find solution X
            XthisNew=SinvAtb-SumInv*Lambda;
            XNew=X;
            XNew(Idxthis)=XthisNew;
            X=XNew;
            break;
        end
    else
        ANewthis=A(:,ActiveNew);
        XNewthis=XNew(ActiveNew,:);
        if LNew<=0 && norm(ANewthis*XNewthis-b,2)>Lambda
            X=repmat(nan,size(X,1),size(X,2));
            fprintf('\nNo Answer!\n');
            return;
        end
        if norm(ANewthis*XNewthis-b,2)<Lambda
            Athis=A(:,Active);
            Xthis=X(Active,:);
            Idxthis=find(Active);
            Inv=(Athis'*Athis)^-1;
            % calculate DLambda and LambdaFinal
            Alpha=Athis*Xthis-b;
            Beta=Athis*Inv*c(Active,:);
            Dot=Alpha'*Beta;
            LAlpha=Alpha'*Alpha;
            LBeta=Beta'*Beta;
            DLambda=(-2*Dot-sqrt(max(0,4*Dot*Dot-4*LBeta*(LAlpha-Lambda*Lambda)))) /(2*LBeta);
            Lfinal=L-DLambda;
            % find the solution X
            SumInv=Inv*c(Active,:);
            Atb=Athis'*b;
            SinvAtb=Inv*Atb;
            XthisNew=SinvAtb-SumInv*Lfinal;
            XNew=X;
            XNew(Idxthis)=XthisNew;
            X=XNew;
            break;
        end
    end
    X=XNew;
    L=LNew;
    Active=ActiveNew;
    [XNew,ActiveNew,LNew]=LassoIteration(A,b,c,L,X,Active,Zero);
    Node.X=XNew;
    Node.L=LNew;
    Node.Active=ActiveNew;
    LassoPath=[LassoPath;cell(1,1)];
    LassoPath{end}=Node;
    LNew

end




end


function [XNew,ActiveNew,LambdaNew]=LassoIteration(A,b,c,Lambda,X,Active,Zero)
nthis=sum(Active>0);
Athis=A(:,Active);
Idxthis=find(Active);
Inv=(Athis'*Athis)^-1;
SumInv=Inv*c(Active,1);
Atb=Athis'*b;
SinvAtb=Inv*Atb;
Xthis=SinvAtb-SumInv*Lambda;
%Check which in Active to become inactive.
LambdaI=-1;
IdxI=0;
for i=1:nthis
    if SumInv(i)>=0
        continue;
    end
    Athis(:,i);
    Lambdathis=Lambda+Xthis(i)/SumInv(i);
    if Lambdathis>LambdaI
        LambdaI=Lambdathis;
        IdxI=Idxthis(i);
    end
end
%Check which out Active to become active.
LambdaA=-1;
IdxA=0;
for i=1:size(A,2)
    if Active(i)
        continue;
    end
    Idxjshort=find(A(:,i)>Zero);
    Ajshort=[A(Idxjshort,i)]';
    AjtA=Ajshort*Athis(Idxjshort,:);
    D=AjtA*SumInv-c(i);
    if D>=0
        continue;
    end
    N=AjtA*Xthis-[A(:,i)]'*b+c(i)*Lambda;
    Lambdathis=Lambda+N/D;
    if Lambdathis>LambdaA
        LambdaA=Lambdathis;
        IdxA=i;
    end
end

ActiveNew=Active;
if LambdaI>LambdaA
% remove 1 element
    if LambdaI>=0
        LambdaNew=LambdaI;
    else
        LambdaNew=0;
    end
    XthisNew=SinvAtb-SumInv*LambdaNew;
    XNew=X;
    XNew(Idxthis)=XthisNew;
    if IdxI~=0
        ActiveNew(IdxI)=0;
        XNew(IdxI)=0;
    end

else
    if LambdaA>=0
        LambdaNew=LambdaA;
    else
        LambdaNew=0;
    end
    if IdxA~=0
        ActiveNew(IdxA)=1;
    end
    XthisNew=SinvAtb-SumInv*LambdaNew;
    XNew=X;
    XNew(Idxthis)=XthisNew;
end


end











