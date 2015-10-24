#include"EMfunctions.h"

void Sparse_malloc(struct Sparse *Sp, int n)
{
    Sp->n=n;
    Sp->X=(short *)malloc(sizeof(short)*n);
    Sp->Y=(short *)malloc(sizeof(short)*n);
    Sp->V=(double *)malloc(sizeof(double)*n);
    return;
}

void Sparse_free(struct Sparse *Sp)
{
    free(Sp->X);
    free(Sp->Y);
    free(Sp->V);
    return;
}

void State_malloc(struct State *St, int n, int NumNo, bool EvenNoise)
{
    St->n=n;
    St->X=(double *)malloc(sizeof(double)*n);
    St->Y=(double *)malloc(sizeof(double)*n);
    St->I=(double *)malloc(sizeof(double)*n);
    if(!EvenNoise){
	St->NumNo=NumNo;
	St->No=(double *)malloc(sizeof(double)*NumNo);
    }
    return;
}

void State_free(struct State *St, bool EvenNoise)
{
    free(St->X);
    free(St->Y);
    free(St->I);
    if(!EvenNoise){
	free(St->No);
    }
    return;
}

void RUN_free(struct RUN *run, bool EvenNoise)
{
    int i;
    free(run->bb);
    free(run->bound);
    free(run->Wn);
    free(run->S);
    free(run->exptable);
    free(run->expx);
    free(run->X);
    free(run->Y);
    free(run->die);
    free(run->freeze);
    for(i=0;i<run->n;i++){
	free(run->W[i].X);
	free(run->W[i].Y);
	free(run->W[i].V);
	free(run->w[i].X);
	free(run->w[i].Y);
	free(run->w[i].V);
    }
    free(run->W);
    free(run->w);
    if(!EvenNoise){
	for(i=0;i<run->NumNo;i++){
	    free(run->NoDist[i].X);
	    free(run->NoDist[i].Y);
	    free(run->NoDist[i].V);
	    free(run->NoAssign[i].X);
	    free(run->NoAssign[i].Y);
	    free(run->NoAssign[i].V);
	}
	free(run->NoDist);
	free(run->NoAssign);
    }
    return;
}

void CPtable_malloc(struct CPtable *CP, int NumLambda, int NumMu)
{
    CP->NumMu=NumMu;
    CP->NumLambda=NumLambda;
    CP->MuX=(double *)malloc(sizeof(double)*NumMu);
    CP->LambdaX=(double *)malloc(sizeof(double)*NumLambda);
    CP->LogMuX=(double *)malloc(sizeof(double)*NumMu);
    CP->LogLambdaX=(double *)malloc(sizeof(double)*NumLambda);
    int i;
    CP->cptable=(double **)malloc(sizeof(double *)*NumLambda);
    for(i=0;i<NumLambda;i++){
        CP->cptable[i]=(double *)malloc(sizeof(double)*NumMu);
    }
    return;
}

void CPtable_free(struct CPtable *CP)
{
    free(CP->MuX);
    free(CP->LambdaX);
    free(CP->LogMuX);
    free(CP->LogLambdaX);
    int i;
    for(i=0;i<CP->NumLambda;i++){
        free(CP->cptable[i]);
    }
    free(CP->cptable);
    return;
}

void CPtable_prepare(struct CPtable *CP)
{
    int i;
    for(i=0;i<CP->NumLambda;i++){
	CP->LogLambdaX[i]=log(CP->LambdaX[i]);
    }
    for(i=0;i<CP->NumMu;i++){
	CP->LogMuX[i]=log(CP->MuX[i]);
    }
    CP->LogMinLambda=CP->LogLambdaX[0];
    CP->LogMaxLambda=CP->LogLambdaX[CP->NumLambda-1];
    CP->LogStepLambda = (CP->LogMaxLambda-CP->LogMinLambda) / (CP->NumLambda-1.0);
    CP->LogMinMu=CP->LogMuX[0];
    CP->LogMaxMu=CP->LogMuX[CP->NumMu-1];
    CP->LogStepMu = (CP->LogMaxMu-CP->LogMinMu) / (CP->NumMu-1.0);
    return;
}

/* function value searched by linear interpolation */
double  calibrated_poisson(double Lambda, double Mu, struct CPtable *CP)
{
    double LogMu,LogLambda;
    LogMu=log(Mu);
    if(LogMu < CP->LogMinMu){
	LogMu=CP->LogMinMu;
    }
    else if(LogMu > CP->LogMaxMu){
	LogMu = CP->LogMaxMu - DoubleZero;
    }
    LogLambda=log(Lambda);
    if(LogLambda < CP->LogMinLambda){
	LogLambda=CP->LogMinLambda;
    }
    else if(LogLambda > CP->LogMaxLambda){
	LogLambda = CP->LogMaxLambda - DoubleZero;
    }
    int MuCell, LambdaCell;
    MuCell=(int)((LogMu-CP->LogMinMu)/CP->LogStepMu);
    LambdaCell=(int)((LogLambda-CP->LogMinLambda)/CP->LogStepLambda);
    double LogMuInCell,LogLambdaInCell;
    LogMuInCell = (LogMu - CP->LogMuX[MuCell]) / CP->LogStepMu;
    LogLambdaInCell = (LogLambda - CP->LogLambdaX[LambdaCell]) / CP->LogStepLambda;
    double ans;
    ans= CP->cptable[LambdaCell][MuCell]*(1-LogLambdaInCell)*(1-LogMuInCell) + CP->cptable[LambdaCell][MuCell+1]*(1-LogLambdaInCell)*LogMuInCell
	+ CP->cptable[LambdaCell+1][MuCell]*LogLambdaInCell*(1-LogMuInCell) + CP->cptable[LambdaCell+1][MuCell+1]*LogLambdaInCell*LogMuInCell;
    return ans;
}


/*Here A should be pre-allocated*/
void PSFGauss(struct Sparse *A, struct RUN *run, int L1, int U1, int L2, int U2, double x, double y, double sigma, int bdecay) 
{
    int i,j,k;
    int cx,cy,l1,l2,u1,u2;
    cx=(int)(floor(x)+0.5);
    cy=(int)(floor(y)+0.5);
    l1=MAX(L1,cx-bdecay);
    l2=MAX(L2,cy-bdecay);
    u1=MIN(U1+1,cx+bdecay+2);
    u2=MIN(U2+1,cy+bdecay+2);
    A->n=0;
    if(u1-l1<1){
	return;
    }
    if(u2-l2<1){
	return;
    }
    double a,b;
    double Dx,Dy;
    double dxinv=(double)Expinterval;
    a=1.0/(2*M_PI*sigma*sigma);
    b=1.0/(2*sigma*sigma);
    int D;
    int expup=Expmin*Expinterval-1;
    double d,dd;
    for(i=l1;i<u1;i++){
	for(j=l2;j<u2;j++){
	    A->X[A->n]=i;
	    A->Y[A->n]=j;
	    Dx=i-x;
	    Dy=j-y;
	    d=-b*(Dx*Dx+Dy*Dy);
	    D=(int)(-d*dxinv);
	    if(D>expup){
		A->V[A->n]=0.0;
	    }
	    else{
		dd=d-run->expx[D];
		A->V[A->n]=a*run->exptable[D]*(1.0+dd);
	    }
	    (A->n)++;
	}
    }
    return;
}


/*no is the intensity of noise, scale it here to get better speed*/
void PSFno(double *An, int sbs, double no)
{
    double x=no/sbs;
    sbs-=1;
    for(;sbs>-1;sbs-=1)
        An[sbs]=x;
    return;
}

/*Compute noise PSF for uneven noise distribution*/
/*NoDist are the noise particles, No are the intensities of noise particles*/
/*An, NoAssign are pre-allocated*/
void UnevenPSFno(double *An, struct Sparse *NoAssign, int sbs, int sb1, int sb2, int NumNo, struct Sparse *NoDist, double *No)
{
    int i,j;
    int idx;
    for(i=0;i<sbs;i++)
        An[i]=0;
    for(i=0;i<NumNo;i++){
	NoAssign[i].n=NoDist[i].n;
	for(j=0;j<NoDist[i].n;j++){
	    idx=NoDist[i].X[j]+NoDist[i].Y[j]*sb1;
	    NoAssign[i].X[j]=NoDist[i].X[j];
	    NoAssign[i].Y[j]=NoDist[i].Y[j];
	    NoAssign[i].V[j]=NoDist[i].V[j]*No[i];
	    An[idx]+=NoAssign[i].V[j];
	}
    }
    double *Aninv;
    Aninv=(double *)malloc(sizeof(double)*sbs);
    for(i=0;i<sbs;i++)
	Aninv[i]=1.0/An[i];
    for(i=0;i<NumNo;i++){
	for(j=0;j<NoDist[i].n;j++){
	    idx=NoDist[i].X[j]+NoDist[i].Y[j]*sb1;
	    NoAssign[i].V[j]=NoAssign[i].V[j]*Aninv[idx];
	}
    }
    free(Aninv);
    return;
}

/*Run one EM iteration.*/ 
/*If PositionFix==true, then Molecule positions will not upgrade*/
void RunStep(struct RUN *run, struct State *pic, struct State *pic0, bool PositionFix, bool DoCalibration, struct CPtable *CP, bool EvenNoise)
{
    int i,j;
    int n=run->n;
    int s1=run->s1;
    int s2=run->s2;
    int bsize=run->bsize;
    double *bb;
    bool *bound=run->bound;
    double *S=run->S;
    struct Sparse *W=run->W;
    double *Wn=run->Wn;
    int sb1=run->sb1;
    int sb2=run->sb2;
    int sbs=run->sbs;
    bool *die=run->die;
    bool *freeze=run->freeze;
/*Initialize*/
    bb=(double *)malloc(sizeof(double)*sbs);
    for(i=0;i<sbs;i++)
	bb[i]=run->bb[i];
/*Estep*/
    /*reset the picture and S*/
    for(i=0;i<sbs;i++){
	S[i]=0;
    }
    double mp;
    struct Sparse *sp,*ssp;
    /*PSF calculation*/
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	if(freeze[i]){
	}
	else{
	    PSFGauss(run->w+i, run,0,sb1-1,0,sb2-1,pic->X[i],pic->Y[i],run->sig,run->psfdecay);
	}
	ssp=run->w+i;
	sp=W+i;
	sp->n=ssp->n;
	for(j=0;j<ssp->n;j++){
	    sp->X[j]=ssp->X[j];
	    sp->Y[j]=ssp->Y[j];
	    sp->V[j]=ssp->V[j]*pic->I[i];
	}
    }
    /*Sum up PSFs*/
    double c;
    int ind;
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	sp=W+i;
	for(j=0;j<sp->n;j++){
	    ind=sp->X[j]+sp->Y[j]*sb1;
	    c=sp->V[j];
	    S[ind]+=c;
	    if(bound[ind]){
		continue;
	    }
	}
    }
    if(EvenNoise){
        PSFno(Wn,sbs,pic->no);
    }
    else{
	UnevenPSFno(Wn,run->NoAssign,sbs,sb1,sb2,pic->NumNo,run->NoDist,pic->No);
    }
    /*reconstruct the hidden part of the picture*/
    double *Sinv;
    Sinv=(double *)malloc(sizeof(double)*sbs);
    double nosum=0;
    for(i=0;i<sbs;i++){
        S[i]+=Wn[i];
	if(bound[i]){
	    bb[i]=S[i];
	}
	else{
	    if(DoCalibration){
		bb[i]=calibrated_poisson(S[i],bb[i],CP);
	    }
	}
	Sinv[i]=1.0/S[i];
	Wn[i]*=Sinv[i]*bb[i];
        nosum+=Wn[i];
    }
    
    /*update W,Wn to image based PSF*/
    double *O,*Lx,*Ly;
    O=(double *)malloc(sizeof(double)*n);
    Lx=(double *)malloc(sizeof(double)*n);
    Ly=(double *)malloc(sizeof(double)*n);
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	sp=W+i;
        O[i]=0;
        Lx[i]=0;
	Ly[i]=0;
	for(j=0;j<sp->n;j++){
	    ind=sp->X[j]+sp->Y[j]*sb1;
	    sp->V[j]*=Sinv[ind]*bb[ind];
            O[i]+=sp->V[j];
	    Lx[i]+=sp->V[j]*sp->X[j];
	    Ly[i]+=sp->V[j]*sp->Y[j];
	}
    }
/*copy pic to pic0*/
    pic0->no=pic->no;
    for(i=0;i<run->n;i++){
	if(die[i])
	    continue;
	pic0->X[i]=pic->X[i];
	pic0->Y[i]=pic->Y[i];
	pic0->I[i]=pic->I[i];
    }
/*Mstep*/
    /*update intensity(np,nosum, O[n])*/
    double Iinv=1.0/(1.0+run->lambda);
    if(EvenNoise){
        pic->no=nosum;
    }
    else{
	for(i=0;i<run->NumNo;i++){
	    pic->No[i]=0;
	}
        for(i=0;i<run->NumNo;i++){
	    for(j=0;j<run->NoAssign[i].n;j++){
		ind=run->NoAssign[i].X[j]+run->NoAssign[i].Y[j]*sb1;
		pic->No[i]+=Wn[ind]*run->NoAssign[i].V[j];
	    }
	}
    }
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	pic->I[i]=O[i]*Iinv;
    }
    /*update locations(Tx[n],Ty[n])*/
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	if(O[i]<Zero){
	    die[i]=true;
	    continue;
	}
	if(PositionFix){
	}
	else
	{
	    pic->X[i]=Lx[i]/O[i];
	    pic->Y[i]=Ly[i]/O[i];
	}
    }
/*Check stop criterion*/
    double errI=fabs(pic->no-pic0->no),errP=0,Norm_errP=0;
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	errI+=fabs(pic->I[i]-pic0->I[i]);
	if(PositionFix){
	}
	else{
	    errP+=pic0->I[i]*(fabs(pic->X[i]-pic0->X[i])+fabs(pic->Y[i]-pic0->Y[i]));
	    Norm_errP+=pic0->I[i];
	}
	mp=(fabs(run->X[i]-pic->X[i])+fabs(run->Y[i]-pic->Y[i]));
	if(mp>MZero){
	    run->X[i]=pic->X[i];
	    run->Y[i]=pic->Y[i];
	    freeze[i]=false;
	}
	else{
	    freeze[i]=true;
	}
    }
    if(PositionFix){
    }
    else{
	errP=errP/Norm_errP;
    }
    if(errI<run->Izero&&errP<run->Pzero){
	run->terminate='Y';
    }
    free(bb);
    free(Sinv);
    free(O);
    free(Lx);
    free(Ly);
    return;
}




