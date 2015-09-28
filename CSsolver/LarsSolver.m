function [ X, LassoPath ] = LarsSolver( A, b, c, Lambda, Type, X0, NoPrintPath )
%Solves this Problem by Least angle regression: 
%    Type=0: 0.5*|(Ax-b)|^2+Lambda*c'*x  s.t. x>=0
%    Type=1: c'*x  s.t. |(Ax-b)|^2 <= Lambda^2, && x>=0

% A,b,c and Lambda are assumed to be non-negative.
%   
% Initial point x0 could be given outside. x0 should be the solution of the
% following optimization problem:
%   minimize: 0.5*|(Ax-b)|^2  s.t. x>=0, x_i=0 if c_i~=0
% If x0 isn't given, it is computed by function quadprog().

Zero=0; % minimal value of A, below which regarded as zero
Czero=1e-15; % minimal value of c
Qzero=1e-10; % value to judge active elements in x0
m=size(A,1);
n=size(A,2);
ActiveNew=zeros(n,1);
ActiveNew=ActiveNew>0.5;
XNew=zeros(n,1);
vpadigit=50

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
    ActiveNew(X0>=Qzero,:)=1;
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
Node.X=sparse(X);
Node.L=L;
Node.Active=sparse(Active);
LassoPath=cell(1,1);
LassoPath{1}=Node;

%%% Run lars iterations
while 1
    if Type==0
        if LNew<Lambda
            Athis=A(:,Active);
            Idxthis=find(Active);
            Atb=Athis'*b;
            
            [U,D,V]=svd(Athis'*Athis);
            if cond(D)==inf || cond(D)>1e12
                fprintf('\nConvert to high pricision algibra\n');
                [U,D,V]=svd(vpa(Athis'*Athis,vpadigit));
                Inv=(U*inv(D)*V');
                SumInv=double(Inv*c(Active,:));
                SinvAtb=double(Inv*Atb);
            else
                Inv=U*inv(D)*V';
                SumInv=Inv*c(Active,:);
                SinvAtb=Inv*Atb;
            end

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
            Atb=Athis'*b;
            
            [U,D,V]=svd(Athis'*Athis);
            if cond(D)==inf || cond(D)>1e12
                fprintf('\nConvert to high pricision algibra\n');
                [U,D,V]=svd(vpa(Athis'*Athis,vpadigit));
                Inv=(U*inv(D)*V');
                Beta=double(Athis*(Inv*c(Active,:)));
                SumInv=double(Inv*c(Active,:));
                SinvAtb=double(Inv*Atb);
            else
                Inv=U*inv(D)*V';
                Beta=Athis*Inv*c(Active,:);
                SumInv=Inv*c(Active,:);
                SinvAtb=Inv*Atb;
            end
            % calculate DLambda and LambdaFinal
            Alpha=Athis*Xthis-b;
            Dot=Alpha'*Beta;
            LAlpha=Alpha'*Alpha;
            LBeta=Beta'*Beta;
            DLambda=(-2*Dot-sqrt(max(0,4*Dot*Dot-4*LBeta*(LAlpha-Lambda*Lambda)))) /(2*LBeta);
            Lfinal=L-DLambda;
            % find the solution X
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
    [XNew,ActiveNew,LNew]=LassoIteration(A,b,c,L,X,Active,Zero,0,vpadigit);
    if LNew>L
        fprintf('\nTerminate before stop condition met\n');
        return;
    end
    Node.X=sparse(XNew);
    Node.L=LNew;
    Node.Active=sparse(ActiveNew);
    LassoPath=[LassoPath;cell(1,1)];
    LassoPath{end}=Node;
    if exist('NoPrintPath')
    else
        LNew
    end

end




end


function [XNew,ActiveNew,LambdaNew]=LassoIteration(A,b,c,Lambda,X,Active,Zero,IfVpa,vpadigit)
nthis=sum(Active>0);
Athis=A(:,Active);
Idxthis=find(Active);
Atb=Athis'*b;
if IfVpa
    [U,D,V]=svd(vpa(Athis'*Athis,vpadigit));
    Inv=(U*inv(D)*V');
    SumInv=double(Inv*c(Active,1));
    SinvAtb=double(Inv*Atb);
else
    [U,D,V]=svd(Athis'*Athis);
    if cond(D)==inf || cond(D)>1e12
        fprintf('\nConvert to high pricision algibra,vpadigit=%d\n',vpadigit);
        [U,D,V]=svd(vpa(Athis'*Athis,vpadigit));
        Inv=(U*inv(D)*V');
        SumInv=double(Inv*c(Active,1));
        SinvAtb=double(Inv*Atb);
    else
        Inv=U*inv(D)*V';
        SumInv=Inv*c(Active,1);
        SinvAtb=Inv*Atb;
    end
end
Xthis=SinvAtb-SumInv*Lambda;
%Check which in Active to become inactive.
IdxMinus=find(Xthis<0);
if length(IdxMinus)==0
else
    DeltaLambda=min(Xthis(IdxMinus,:)./SumInv(IdxMinus,:));
    Lambda=Lambda+DeltaLambda;
    Xthis=SinvAtb-SumInv*Lambda;
end
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
    N=AjtA*Xthis-[A(:,i)]'*b+c(i)*Lambda;
    if D>=0
        continue;
    end
    if N/D>0
        LambdaNew=Lambda+1;
        XNew=X;
        ActiveNew=Active;
        return;
    end
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
    XNew=X;
    if IdxI~=0
        ActiveNew(IdxI)=0;
        XNew(IdxI)=0;
        % calibrate lambda
        AthisNew=A(:,ActiveNew);
        IdxthisNew=find(ActiveNew);
        AtbNew=AthisNew'*b;
        if IfVpa
            [U,D,V]=svd(vpa(AthisNew'*AthisNew,vpadigit));
            InvNew=(U*inv(D)*V');
            SumInvNew=double(InvNew*c(ActiveNew,1));
            SinvAtbNew=double(InvNew*AtbNew);
        else
            [U,D,V]=svd(AthisNew'*AthisNew);
            if cond(D)==inf || cond(D)>1e12
                fprintf('\nConvert to high pricision algibra,vpadigit=%d\n',vpadigit);
                [U,D,V]=svd(vpa(AthisNew'*AthisNew,vpadigit));
                InvNew=(U*inv(D)*V');
                SumInvNew=double(InvNew*c(ActiveNew,1));
                SinvAtbNew=double(InvNew*AtbNew);
            else
                InvNew=U*inv(D)*V';
                SumInvNew=InvNew*c(ActiveNew,1);
                SinvAtbNew=InvNew*AtbNew;
            end
        end
        XthisNew=SinvAtbNew-SumInvNew*LambdaNew;
        IdxjshortNew=find(A(:,IdxI)>Zero);
        AjshortNew=[A(IdxjshortNew,IdxI)]';
        AjtANew=AjshortNew*AthisNew(IdxjshortNew,:);
        D=AjtANew*SumInvNew-c(IdxI);
        N=AjtANew*XthisNew-[A(:,IdxI)]'*b+c(IdxI)*LambdaNew;
        if N<0
            LambdaNew=(AjtANew*SinvAtbNew-[A(:,IdxI)]'*b)/D;
        end
        XthisNew=SinvAtbNew-SumInvNew*LambdaNew;
        XNew(IdxthisNew)=XthisNew;
    else
        XthisNew=SinvAtb-SumInv*LambdaNew;
        XNew(Idxthis)=XthisNew;
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












