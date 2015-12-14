function [ Losses, LossesDetail, BestLambdaIdx, KLlosses ] = LossKL( Img, LambdaPath, RLasso, psfno, OptPara, DoTunning)
%LOSSKL Summary of this function goes here
%   Detailed explanation goes here
L=length(RLasso);
Losses=zeros(L,1);
LossesDetail=zeros(L,2);
Zero=1e-12;
Tolerance=[5000,1e-6*sum(sum(Img)),1e-5];
Mergedist=0.1;
MolZero1=0.5;
MolZero2=10;
MolZero3=200;
FoldRemain=0.9999;
LambdaLocal=0.1;
Nsample=100;
FoldRelax=0;
LambdaOrigin=0.2;
KLlosses=cell(L,1);
for i=1:L
    if DoTunning==1
        mv0.pic=RLasso{i}.pic;
        mv0.no=RLasso{i}.no;
        [Pic,No]=RunSolverFix(Img, OptPara, Tolerance, mv0);
    else
        Pic=RLasso{i}.pic;
        No=RLasso{i}.no;
    end
    Pic=PostRun(Pic,[Mergedist,MolZero1,FoldRemain]);
    Pic=PostRun(Pic,[Mergedist,MolZero2,FoldRemain]);
    [ Signal, Bg, StandardPSF, ExpectedSignal, ExpectedBg ] = SignalSeparate( Img, Pic, No, psfno, OptPara );
    ExpectedFull=ExpectedBg+sum(ExpectedSignal,2);
    KL=zeros(size(Signal,2),1);
    KL1=zeros(size(Signal,2),1);
    KL2=zeros(size(Signal,2),1);
    Weight=Pic(:,3)-(Pic(:,3)>MolZero3).*(Pic(:,3)-MolZero3);
    Weight=Weight/MolZero3;
    for j=1:size(Signal,2)
        Pa=full(ExpectedSignal(:,j));
        Pb=ExpectedFull-Pa;
        P1=full(Signal(:,j));
        P2=full(StandardPSF(:,j));
        P2=P2/sum(P2);
        v1=find(P1>Zero);
        v2=find(P2>Zero);
        E=zeros(1,Nsample);
        P3=zeros(size(P1,1),Nsample);
        Psig=2*sqrt(Pa(v1,:).*Pb(v1,:)./(Pa(v1,:)+Pb(v1,:)));
        P3(v1,:)=repmat(P1(v1,:),1,Nsample)+repmat(Psig,1,Nsample).*randn(length(v1),Nsample);
        P3=(P3+abs(P3))/2;
        P3=P3/sum(P1);
        %P3=P3./repmat(sum(P3,1),size(P3,1),1);
        TA1=zeros(size(P3));
        v3=find(P3>Zero);
        TA1(v3)=P3(v3).*log(P3(v3));
        TA2=P3(v2,:).*log(repmat(P2(v2,:),1,Nsample));
        E=sum(TA1,1)-sum(TA2,1);
        P1=P1/sum(P1);
        KL1(j)=[P1(v1)]'*log(P1(v1))-[P1(v2)]'*log(P2(v2));
        TA3=P3(v1,:).*log(repmat(P1(v1,:),1,Nsample));
        KL2(j)=mean(sum(TA1,1)-sum(TA3,1));
        P1=sqrt(P1);
        P1=P1/sum(P1);
        P2=sqrt(P2);
        P2=P2/sum(P2);
        KL(j)=[P1(v1)]'*log(P1(v1))-[P1(v2)]'*log(P2(v2));
    end
    KLlosses{i}.KL=KL;
    KLlosses{i}.KL1=KL1;
    KLlosses{i}.KL2=KL2;
    Losses(i)=sum(Weight.*KL)/sum(Weight);
    LossesDetail(i,1)=sum(Weight.*KL1)/sum(Weight);
    LossesDetail(i,2)=sum(Weight.*KL2)/sum(Weight);
end

% find best lambda
  % first minimum
for BestLambdaIdx=1+1:L
    ThisIs=true;
    for i=BestLambdaIdx+1:L
        if LambdaPath(i)-LambdaPath(BestLambdaIdx)>LambdaLocal
            break;
        end
        if Losses(i)*(1+FoldRelax)<Losses(BestLambdaIdx);
            ThisIs=false;
            break;
        end
    end
    if ThisIs==true
        break;
    end
    if LambdaPath(L)-LambdaPath(BestLambdaIdx)<LambdaLocal
        BestLambdaIdx=0;
        break;
    end
end

  % second minimum
if BestLambdaIdx==0
    return;
end
LossesRemain=Losses(BestLambdaIdx+1:end);
v=find(LossesRemain<Losses(BestLambdaIdx));
    % no second minimum
if length(v)==0
    return;
end
    % find second minimum
for k=1:length(v)
    BestLambdaIdx2=v(k)+BestLambdaIdx;
    ThisIs=true;
    for i=BestLambdaIdx2+1:L
        if LambdaPath(i)-LambdaPath(BestLambdaIdx2)>LambdaLocal
            break;
        end
        if Losses(i)*(1+FoldRelax)<Losses(BestLambdaIdx2);
            ThisIs=false;
            break;
        end
    end
    if ThisIs==true
        break;
    end
    if LambdaPath(L)-LambdaPath(BestLambdaIdx2)<LambdaLocal
        BestLambdaIdx2=0;
        break;
    end
end
    % judge which minimum is better
if BestLambdaIdx2==0
    return;
else
    if Losses(BestLambdaIdx)*(LambdaOrigin+LambdaPath(BestLambdaIdx))<Losses(BestLambdaIdx2)*(LambdaOrigin+LambdaPath(BestLambdaIdx2))
    else
        BestLambdaIdx=BestLambdaIdx2;
    end
end

end

