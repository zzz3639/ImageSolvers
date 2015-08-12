/* Modified in 2015.07.10 by ZHANG Haowen
   Do SparseGMM fitting by EM algorithm
   Positions are fixed
*/
#include"mex.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define Zero 1e-1
#define MZero 1e-8

#define Expmin 32
#define Expinterval 1024 

struct Sparse{
    int n;
    short *X;
    short *Y;
    double *V;
};

struct InType{
    int s1;
    int s2;
    double *b;
    int n;
    double sig;
    int MaxIt;
    double Pzero;
    double Izero;
    double lambda;
    int bsize;
    int bdecay;
    long int interval;
    long int seed;
};
struct State{
    double no;
    double *X;
    double *Y;
    double *I;
};
struct RUN{
    int n;
    int s1;
    int s2;
    int sb1;
    int sb2;
    int ss;
    int sbs;
    double *bb;  /*Image with a bdecay size boundary*/
    bool *bound; /*Belongs to boundary or not*/
    struct Sparse *W;  /*Record PSF*/
    struct Sparse *w;
    double *Wn;  /*PSF of noise*/
    double *S;  
    double Int;
    char terminate;
    double *exptable; /*Accelerate exp function*/
    double *expx;
    double *X;  /*spots locations from last iteration*/
    double *Y;
    bool *die;  /*spots with an intensity zero*/
    bool *freeze; /*spots didn't move from past few iterations*/
};
struct MV{
    int T;
    int P;
    double **X;
    double **Y;
    double **I;
    double *no;
};




/*Here A should be pre-allocated*/
void PSF(struct Sparse *A, struct RUN *run, int L1, int U1, int L2, int U2, double x, double y, double sigma, int bdecay) 
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

void ReadInput(int nrhs ,const mxArray *prhs[], struct InType *input)
{
    int k;
/*reading b;*/
    int s1,s2;
    s1=mxGetM(prhs[0]);
    s2=mxGetN(prhs[0]);
    input->s1=s1;
    input->s2=s2;
    int i,j;
    double *t;
    t=mxGetPr(prhs[0]);
    input->b=(double *)malloc(sizeof(double)*s1*s2);
    for(i=0;i<s1*s2;i++){
        *(input->b+i)=*(t+i);
    }
/*reading n;*/
    t=mxGetPr(prhs[1]);
    input->n=(int)(*t+0.5);
/*reading sig;*/
    input->sig=*(mxGetPr(prhs[2]));
/*reading stop criterions*/  
    t=mxGetPr(prhs[3]);
    input->MaxIt=(int)(*t+0.5);
    input->Izero=*(t+1);
    input->Pzero=*(t+2);
/*reading boundary size;*/
    input->bsize=(int)(*(mxGetPr(prhs[4]))+0.5);
/*reading PSF range;*/
    input->bdecay=(int)(*(mxGetPr(prhs[5]))+0.5);
/*reading lambda;*/
    input->lambda=*(mxGetPr(prhs[6]));
/*reading interval length*/
    t=mxGetPr(prhs[7]);
    input->interval=*((long int *)t);
    input->seed=-1;
    if(mxGetN(prhs[7])>1)
	input->seed=*((long int *)t+1);
    if(nrhs>8){
    }
    return;
}

void WriteOutput(int nlhs, mxArray *plhs[], int n, struct State *pic, struct MV *mv, long int seed, int bsize)
{
    if(nlhs<1)
        return;
/*Write pic and no to the output*/
    plhs[0]=mxCreateNumericMatrix(n,3,mxDOUBLE_CLASS,mxREAL);
    double *point;
    point=mxGetPr(plhs[0]);
    int i,j;
    for(i=0;i<n;i++){
	*(point+i)=pic->Y[i]-bsize;
	*(point+i+n)=pic->X[i]-bsize;
	*(point+i+n+n)=pic->I[i];
    }
    if(nlhs<2)
	return;
    plhs[1]=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(plhs[1]);
    *(point)=pic->no;

    if(nlhs<3)
        return;
/*Write MV to the output*/
    const char *fn1="T";
    const char **fnames=&fn1;
    plhs[2]=mxCreateStructMatrix(1, 1, 1, fnames);
    mxArray *pv;
    pv=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(pv);
    *point=mv->T;
    mxSetField(plhs[2],0,"T",pv);

    mxAddField(plhs[2],"N");
    pv=mxCreateNumericMatrix(1,1,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(pv);
    *point=mv->P;
    mxSetField(plhs[2],0,"N",pv);

    mxAddField(plhs[2],"no");
    pv=mxCreateNumericMatrix(mv->P,1,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(pv);
    for(i=0;i<mv->P;i++)
        *(point+i)=mv->no[i];
    mxSetField(plhs[2],0,"no",pv);

    mxAddField(plhs[2],"pic");
    const int dms[3]={n,3,mv->P};
    pv=mxCreateNumericArray(3,dms,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(pv);
    for(i=0;i<mv->P;i++){
	for(j=0;j<n;j++){
	    *(point+i*3*n+j)=mv->Y[i][j]-bsize;
	    *(point+i*3*n+j+n)=mv->X[i][j]-bsize;
	    *(point+i*3*n+j+2*n)=mv->I[i][j];
	}
    }
    mxSetField(plhs[2],0,"pic",pv);

    if(nlhs<4)
	return;
    plhs[3]=mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
    point=mxGetPr(plhs[3]);
    *((long int *)point)=seed;

}

void PrintInput(struct InType *input)
{
    mexPrintf("\nInput data and parameters:\n");
    mexPrintf("Picture parameters: s=%d*%d\n",input->s1,input->s2);
    mexPrintf("Optimization parameters: n=%d, sigma=%f, boundary size=%d, PSF range=%d, lambda=%f\n", input->n,input->sig,input->bsize,input->bdecay,input->lambda);
    mexPrintf("Stop criterion: Max iteration number=%d, Intensity stop criterion=%e, Position stop criterion=%e\n",input->MaxIt,input->Izero,input->Pzero);
    mexPrintf("Save every %ld iterations\n",input->interval);
    mexPrintf("Random Seed=%ld\n\n",input->seed);
    return;
}

void Initialize(struct InType *input, struct RUN *R)
{
    mexPrintf("\nStart to initialize\n");
    R->n=input->n;
    R->s1=input->s1;
    R->s2=input->s2;
    R->ss=input->s1*input->s2;
    R->sb1=input->s1+2*input->bsize;
    R->sb2=input->s2+2*input->bsize;
    R->sbs=R->sb1*R->sb2;
    R->bb=(double *)malloc(sizeof(double)*R->sbs);
    R->bound=(bool *)malloc(sizeof(bool)*R->sbs);
    int i,j,a,b;
    for(i=0;i<R->sbs;i++){
	R->bb[i]=0;
	R->bound[i]=true;
    }
    
    int aa;
    R->Int=0;
    for(a=0;a<input->s1;a++){
	aa=a+input->bsize;
        for(b=0;b<input->s2;b++){
	    R->Int+=input->b[a+b*R->s1];
	    R->bb[aa+(input->bsize+b)*R->sb1]=input->b[a+b*R->s1];
	    R->bound[aa+(input->bsize+b)*R->sb1]=false;
	}
    }
    R->S =(double*)malloc(sizeof(double)*R->sbs);
    R->Wn=(double*)malloc(sizeof(double)*R->sbs);
    R->W =(struct Sparse *)malloc(sizeof(struct Sparse)*R->n);
    R->w =(struct Sparse *)malloc(sizeof(struct Sparse)*R->n);
    for(i=0;i<R->n;i++){
	R->W[i].n=0;
	R->W[i].X=(short *)malloc(sizeof(short)*(2+2*input->bdecay)*(2+2*input->bdecay));
	R->W[i].Y=(short *)malloc(sizeof(short)*(2+2*input->bdecay)*(2+2*input->bdecay));
	R->W[i].V=(double *)malloc(sizeof(double)*(2+2*input->bdecay)*(2+2*input->bdecay));
	R->w[i].n=0;
	R->w[i].X=(short *)malloc(sizeof(short)*(2+2*input->bdecay)*(2+2*input->bdecay));
	R->w[i].Y=(short *)malloc(sizeof(short)*(2+2*input->bdecay)*(2+2*input->bdecay));
	R->w[i].V=(double *)malloc(sizeof(double)*(2+2*input->bdecay)*(2+2*input->bdecay));
    }
    R->terminate='N';
    R->exptable=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    R->expx=(double *)malloc(sizeof(double)*Expmin*Expinterval);
    int k;
    double dx;
    k=Expmin*Expinterval;
    dx=1.0/Expinterval;
    for(i=0;i<k;i++){
	R->expx[i]=-dx*i;
	R->exptable[i]=exp(-dx*i);
    }
    mexPrintf("\nInitialize finished\n");

    R->X=(double *)malloc(sizeof(double)*input->n);
    R->Y=(double *)malloc(sizeof(double)*input->n);
    R->die=(bool *)malloc(sizeof(bool)*input->n);
    R->freeze=(bool *)malloc(sizeof(bool)*input->n);
    for(i=0;i<input->n;i++){
	R->X[i]=-1;
	R->Y[i]=-1;
	R->die[i]=false;
	R->freeze[i]=false;
    }
    return;
}

void RunStep(struct RUN *run, struct InType *input, struct State *pic, struct State *pic0)
{
    int i,j;
    int n=input->n;
    int s1=input->s1;
    int s2=input->s2;
    int bsize=input->bsize;
    double *bb=run->bb;
    bool *bound=run->bound;
    double *S=run->S;
    struct Sparse *W=run->W;
    double *Wn=run->Wn;
    int sb1=run->sb1;
    int sb2=run->sb2;
    int sbs=run->sbs;
    bool *die=run->die;
    bool *freeze=run->freeze;
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
	    PSF(run->w+i, run,0,sb1-1,0,sb2-1,pic->X[i],pic->Y[i],input->sig,input->bdecay);
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
    double scale=0;
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
	    scale+=c;
	}
    }
    PSFno(Wn,sbs,pic->no);
    scale+=pic->no/sbs*s1*s2;
    /*reconstruct the hidden part of the picture*/
    scale=run->Int/scale;
    double Sinv[sbs];
    double nosum=0;
    for(i=0;i<sbs;i++){
        S[i]+=Wn[i];
	if(bound[i]){
	    bb[i]=S[i];
	}
	Sinv[i]=1.0/S[i];
	Wn[i]*=Sinv[i]*bb[i];
        nosum+=Wn[i];
    }
    
    /*update W,Wn to image based PSF*/
    double O[n];
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	sp=W+i;
        O[i]=0;
	for(j=0;j<sp->n;j++){
	    ind=sp->X[j]+sp->Y[j]*sb1;
	    sp->V[j]*=Sinv[ind]*bb[ind];
            O[i]+=sp->V[j];
	}
    }
/*copy pic to pic0*/
    pic0->no=pic->no;
    for(i=0;i<input->n;i++){
	if(die[i])
	    continue;
	pic0->I[i]=pic->I[i];
    }
/*Mstep*/
    /*update intensity(np,nosum, O[n])*/
    double Iinv=1.0/(1.0+input->lambda);
    pic->no=nosum;
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	pic->I[i]=O[i]*Iinv;
    }
/*Check stop criterion*/
    double errI=fabs(pic->no-pic0->no);
    for(i=0;i<n;i++){
	if(die[i])
	    continue;
	errI+=fabs(pic->I[i]-pic0->I[i]);
	freeze[i]=true;
    }
    if(errI<input->Izero){
	run->terminate='Y';
    }
    return;
}

void FreeMem(struct RUN *run, struct InType *input, struct MV *mv, struct State *pic, struct State *pic0)
{
/*free mv and pic, pic0;*/
    free(pic0->X);
    free(pic0->Y);
    free(pic0->I);
    free(pic->X);
    free(pic->Y);
    free(pic->I);
    int i;
    for(i=0;i<mv->P;i++){
	free(mv->X[i]);
	free(mv->Y[i]);
	free(mv->I[i]);
    }
    free(mv->X);
    free(mv->Y);
    free(mv->I);
    free(mv->no);
/*free input*/
    free(input->b);
/*free run*/
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
    for(i=0;i<input->n;i++){
	free(run->W[i].X);
	free(run->W[i].Y);
	free(run->W[i].V);
	free(run->w[i].X);
	free(run->w[i].Y);
	free(run->w[i].V);
    }
    free(run->W);
    free(run->w);

    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=9){
        mexErrMsgTxt("\nparameter number incorrect!\nUsage:\n  [pic,no,mv,seed]=EMboundarysparse(b,n,sigma,[MaxIte,Izero,Pzero],bsize,bdecay,Lambda,[int64(NIte)],MV0)\n");
        return;
    }
    struct InType input;
    ReadInput(nrhs,prhs,&input);
    struct RUN run;
    Initialize(&input,&run);

/*Initialize mv*/
    struct MV mv;
    int Smax=(input.MaxIt-1)/input.interval+1;
    mv.X =(double **)malloc(sizeof(double *)*Smax);
    mv.Y =(double **)malloc(sizeof(double *)*Smax);
    mv.I =(double **)malloc(sizeof(double *)*Smax);
    mv.no=(double *)malloc(sizeof(double)*Smax);

/*allocate pic0, pic*/
    struct State pic0,pic;
    pic.X=(double*)malloc(sizeof(double)*input.n);
    pic.Y=(double*)malloc(sizeof(double)*input.n);
    pic.I=(double*)malloc(sizeof(double)*input.n);
    pic0.X=(double*)malloc(sizeof(double)*input.n);
    pic0.Y=(double*)malloc(sizeof(double)*input.n);
    pic0.I=(double*)malloc(sizeof(double)*input.n);

    int i,j,k,t;
/*Initialize pic*/
    {
	mxArray *pv;
	double *mvint;
	pv=mxGetField(prhs[8],0,"no");
	mvint=mxGetPr(pv);
	pic.no=*mvint;
	pv=mxGetField(prhs[8],0,"pic");
	mvint=mxGetPr(pv);
	for(i=0;i<input.n;i++){
	    pic.Y[i]=*(mvint+i)+input.bsize;
	    pic.X[i]=*(mvint+i+input.n)+input.bsize;
	    pic.I[i]=*(mvint+i+input.n+input.n);
	}
    }
    PrintInput(&input);

/*Run EM*/
    mv.P=0;
    mexPrintf("\nRunning\n");
    for(t=0;t<input.MaxIt;t++){
    /*record current snapshot to mv*/
        if(t%input.interval==0){
	    mv.X[mv.P]=(double *)malloc(sizeof(double)*input.n);
	    mv.Y[mv.P]=(double *)malloc(sizeof(double)*input.n);
	    mv.I[mv.P]=(double *)malloc(sizeof(double)*input.n);
	    mv.no[mv.P]=pic.no;
	    for(i=0;i<input.n;i++){
		mv.X[mv.P][i]=pic.X[i];
		mv.Y[mv.P][i]=pic.Y[i];
		mv.I[mv.P][i]=pic.I[i];
	    }
	    mv.P++;
	}
    /*Run EM*/
	RunStep(&run,&input,&pic,&pic0);
    /*Judge break or not*/
        if(run.terminate=='Y')
	    break;
    }
    if(t==input.MaxIt)
        mexPrintf("\nReached maximal iteration number: %d iterations\n\n",t);
    else
	mexPrintf("\nFinished, converge after %d iterations\n\n",t);
    mv.T=t;
    WriteOutput(nlhs,plhs,input.n,&pic,&mv,input.seed,input.bsize);
    FreeMem(&run,&input,&mv,&pic,&pic0);

    return;
}








