function [ Xans ] = ImgTrace( A, b, c, Lambda, Type, LassoPath )
%IMGTRACE Summary of this function goes here
%   X=ImgTrace( A, b, c, Lambda, Type, LassoPath );

vpadigit=50;
L=zeros(length(LassoPath),1);

for i=1:length(LassoPath)
    L(i)=LassoPath{i}.L;
end

%%% track the solution
if Type==0
    % find lambda region
    if Lambda<=L(end)
        Xans=-ones(length(c),1);
        return;
    end
    if Lambda>=L(1)
        Xans=zeros(length(c),1);
        return;
    end
    K=sum(Lambda<L);
    Active=LassoPath{K}.Active;
    X=LassoPath{K}.X;
    % calculate answer
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
    Xans=XNew;
    return;
else
    E=zeros(length(LassoPath),1);
    for i=1:length(LassoPath)
        Active=LassoPath{i}.Active;
        X=LassoPath{i}.X;
        Athis=A(:,Active);
        Xthis=X(Active,:);
        E(i)=norm(Athis*Xthis-b,2);
    end
    if Lambda>=E(1)
        Xans=zeros(length(c),1);
        return;
    end
    if Lambda<=E(end)
        Xans=-ones(length(c),1);
        return;
    end
    K=sum(Lambda<E);
    Active=LassoPath{K}.Active;
    X=LassoPath{K}.X;
    
    Athis=A(:,Active);
    Xthis=X(Active,:);
    Idxthis=find(Active);
    Atb=Athis'*b;
    [U,D,V]=svd(Athis'*Athis);
    if cond(D)==inf || cond(D)>1e12
        fprintf('\nConvert to high pricision algibra\n');
        [U,D,V]=svd(vpa(Athis'*Athis,vpadigit));
        Inv=(U*inv(D)*V');
        Beta=double(Athis*Inv*c(Active,:));
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
    Xans=XNew;
end

end

