/* Modified in 2015.07.10 by ZHANG Haowen
   Do SparseGMM fitting by EM algorithm
   Background noise is modeled by special particles, psf equals to first order B-spline bases.
*/
#include"mex.h"
#include"EMfunctions.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>

#ifndef PosFix
    #define PosFix false
#endif

#ifndef EveNoise
    #define EveNoise false
#endif

struct InType{
    int s1;
    int s2;
    int g1;
    int g2;
    double *b;
    int n;
    int NumNo;
    double sig;
    int MaxIt;
    double Pzero;
    double Izero;
    double lambda;
    int bsize;
    int psfdecay;
    long int interval;
    long int seed;
};

struct MV{
    int T;
    int P;
    double **X;
    double **Y;
    double **I;
    double **No;
};


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
    input->psfdecay=(int)(*(mxGetPr(prhs[5]))+0.5);
/*reading noise grid parameter g1,g2*/
    input->g1=(int)(*(mxGetPr(prhs[6]))+0.5);
    input->g2=(int)(*(mxGetPr(prhs[6])+1)+0.5);
/*reading lambda;*/
    input->lambda=*(mxGetPr(prhs[7]));
/*reading interval length*/
    t=mxGetPr(prhs[8]);
    input->interval=*((long int *)t);
    input->seed=-1;
    if(mxGetN(prhs[8])>1)
	input->seed=*((long int *)t+1);
    if(nrhs>9){
    }
    return;
}

void WriteOutput(int nlhs, mxArray *plhs[], struct RUN *run, struct State *pic, struct MV *mv, long int seed, int bsize)
{
    if(nlhs<1)
        return;
/*Write pic and no to the output*/
    int n,NumNo;
    n=run->n; NumNo=run->NumNo;
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
    plhs[1]=mxCreateNumericMatrix(NumNo,1,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(plhs[1]);
    for(i=0;i<NumNo;i++){
        *(point+i)=pic->No[i];
    }

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
    pv=mxCreateNumericMatrix(mv->P,NumNo,mxDOUBLE_CLASS,mxREAL);
    point=mxGetPr(pv);
    for(i=0;i<mv->P;i++)
	for(j=0;j<NumNo;j++){
	    *(point+i+j*mv->P)=mv->No[i][j];
	}
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
    plhs[3]=mxCreateCellMatrix(NumNo, 1);
    for(i=0;i<NumNo;i++){
	pv=mxCreateNumericMatrix(run->NoDist[i].n,3,mxDOUBLE_CLASS,mxREAL);
	point=mxGetPr(pv);
	for(j=0;j<run->NoDist[i].n;j++){
	    *(point+j)=run->NoDist[i].X[j];
	    *(point+run->NoDist[i].n+j)=run->NoDist[i].Y[j];
	    *(point+2*run->NoDist[i].n+j)=run->NoDist[i].V[j];
	}
	mxSetCell(plhs[3],i,pv);
    }

    if(nlhs<5)
	return;
    plhs[4]=mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
    point=mxGetPr(plhs[4]);
    *((long int *)point)=seed;

}

void PrintInput(struct InType *input)
{
    mexPrintf("\nInput data and parameters:\n");
    mexPrintf("Picture parameters: s=%d*%d\n",input->s1,input->s2);
    mexPrintf("Background noise grid size: g=%d*%d, Noise Particle number: %d\n",input->g1,input->g2, input->NumNo);
    mexPrintf("Optimization parameters: n=%d, sigma=%f, boundary size=%d, PSF range=%d, lambda=%f\n", input->n,input->sig,input->bsize,input->psfdecay,input->lambda);
    mexPrintf("Stop criterion: Max iteration number=%d, Intensity stop criterion=%e, Position stop criterion=%e\n",input->MaxIt,input->Izero,input->Pzero);
    mexPrintf("Save every %ld iterations\n",input->interval);
    mexPrintf("Random Seed=%ld\n\n",input->seed);
    return;
}

void Initialize(struct InType *input, struct RUN *R)
{
    mexPrintf("\nStart to initialize\n");
    R->n=input->n;
    R->lambda=input->lambda;
    R->sig=input->sig;
    R->bsize=input->bsize;
    R->psfdecay=input->psfdecay;
    R->Izero=input->Izero;
    R->Pzero=input->Pzero;
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
	R->W[i].X=(short *)malloc(sizeof(short)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
	R->W[i].Y=(short *)malloc(sizeof(short)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
	R->W[i].V=(double *)malloc(sizeof(double)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
	R->w[i].n=0;
	R->w[i].X=(short *)malloc(sizeof(short)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
	R->w[i].Y=(short *)malloc(sizeof(short)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
	R->w[i].V=(double *)malloc(sizeof(double)*(2+2*input->psfdecay)*(2+2*input->psfdecay));
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

    /*Allocate space to record frozen particles*/
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

    /*Initialize noise particles*/
    int gs1,gs2;
    gs1=(input->g1/2)+R->bsize;
    gs2=(input->g2/2)+R->bsize;
    int numg1,numg2;
    numg1=(R->s1-1-input->g1/2)/input->g1+1;
    numg2=(R->s2-1-input->g2/2)/input->g2+1;

    int *gidx1,*gidx2;
    gidx1=(int *)malloc(sizeof(int)*(numg1+2));
    gidx2=(int *)malloc(sizeof(int)*(numg2+2));
    gidx1[0]=-1;
    for(i=1;i<numg1+1;i++)
	gidx1[i]=gs1+(i-1)*input->g1;
    gidx1[numg1+1]=2*R->bsize+R->s1;
    gidx2[0]=-1;
    for(i=1;i<numg2+1;i++)
	gidx2[i]=gs2+(i-1)*input->g2;
    gidx2[numg2+1]=2*R->bsize+R->s2;

    R->NumNo=numg1*numg2;
    R->NoDist=(struct Sparse *)malloc(sizeof(struct Sparse)*R->NumNo);
    R->NoAssign=(struct Sparse *)malloc(sizeof(struct Sparse)*R->NumNo);
    int LengthPsfnoThis;
    int ig,jg;
    int idxno;
    double mg1,mg2;
    double SumPsfnoThis;
    for(i=0;i<numg1;i++){
	for(j=0;j<numg2;j++){
	    idxno=i+j*numg1;
	    LengthPsfnoThis=(gidx1[i+2]-gidx1[i]-1)*(gidx2[j+2]-gidx2[j]-1);
	    Sparse_malloc(R->NoDist+idxno,LengthPsfnoThis);
	    Sparse_malloc(R->NoAssign+idxno,LengthPsfnoThis);
	    k=0;
	    SumPsfnoThis=0;
	    for(ig=gidx1[i]+1;ig<gidx1[i+2];ig++){
		for(jg=gidx2[j]+1;jg<gidx2[j+2];jg++){
                    R->NoDist[idxno].X[k]=ig;
                    R->NoDist[idxno].Y[k]=jg;
                    R->NoAssign[idxno].X[k]=ig;
                    R->NoAssign[idxno].Y[k]=jg;
		    if(i==0&&ig<gidx1[i+1])
			mg1=1.0;
		    else if(i==numg1-1&&ig>gidx1[i+1])
			mg1=1.0;
		    else
			mg1=1.0-1.0*abs(ig-gidx1[i+1])/input->g1;
		    if(j==0&&jg<gidx2[j+1])
			mg2=1.0;
		    else if(j==numg2-1&&jg>gidx2[j+1])
			mg2=1.0;
		    else
			mg2=1.0-1.0*abs(jg-gidx2[j+1])/input->g2;
		    SumPsfnoThis+=mg1*mg2;
		    R->NoDist[idxno].V[k]=mg1*mg2;
		    k+=1;
		}
	    }
	    for(k=0;k<R->NoDist[idxno].n;k++){
		R->NoDist[idxno].V[k]=R->NoDist[idxno].V[k]/SumPsfnoThis;
	    }
	}
    }

    mexPrintf("\nInitialize finished\n");
    return;
}


void FreeMem(struct RUN *run, struct InType *input, struct MV *mv, struct State *pic, struct State *pic0)
{
/*free mv and pic, pic0;*/
    State_free(pic0,EveNoise);
    State_free(pic,EveNoise);
    int i;
    for(i=0;i<mv->P;i++){
	free(mv->X[i]);
	free(mv->Y[i]);
	free(mv->I[i]);
	free(mv->No[i]);
    }
    free(mv->X);
    free(mv->Y);
    free(mv->I);
    free(mv->No);
/*free input*/
    free(input->b);
/*free run*/
    RUN_free(run,EveNoise);

    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=9&&nrhs!=10){
        mexErrMsgTxt("\nparameter number incorrect!\nUsage:\n [pic,no,mv,PSFNo,seed]=EMUneven(b,n,sigma,[MaxIte,Izero,Pzero],bsize,psfdecay,[g1,g2],Lambda,[int64(NIte)])\nOr\n [pic,no,mv,PSFNo,seed]=EMUneven(b,n,sigma,[MaxIte,Izero,Pzero],bsize,psfdecay,[g1,g2],Lambda,[NIte,Seed])\nOr\n [pic,no,mv,PSFNo,seed]=EMUneven(b,n,sigma,[MaxIte,Izero,Pzero],bsize,psfdecay,[g1,g2],Lambda,[int64(NIte)],MV0)\n");
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
    mv.No=(double **)malloc(sizeof(double *)*Smax);

/*allocate pic0, pic*/
    struct State pic0,pic;
    State_malloc(&pic, input.n, run.NumNo, EveNoise);
    State_malloc(&pic0,input.n, run.NumNo, EveNoise);

    int i,j,k,t;
/*Initialize pic*/
    input.NumNo=run.NumNo;
    {
        if(nrhs==9){
	    if(input.seed<0){
		srand(clock());
		input.seed=time(NULL)+rand();
		srand(input.seed);
	    }
	    else
		srand(input.seed);
	    double IntImgTotal=0;
	    for(i=0;i<input.s1*input.s2;i++){
		IntImgTotal+=input.b[i];
	    }
	    IntImgTotal=IntImgTotal*run.sbs/run.ss;
	    double X1=IntImgTotal*0.5/(input.n);
	    double X2=IntImgTotal*0.5/(run.NumNo);
	    for(i=0;i<input.NumNo;i++){
		pic.No[i]=X2;
	    }
	    for(i=0;i<input.n;i++){
		pic.X[i]=1.0*rand()/RAND_MAX*(input.s1-1)+input.bsize;
		pic.Y[i]=1.0*rand()/RAND_MAX*(input.s2-1)+input.bsize;
		pic.I[i]=X1;
	    }
	}
	else{
	    mxArray *pv;
	    double *mvint;
	    pv=mxGetField(prhs[9],0,"no");
	    mvint=mxGetPr(pv);
	    for(i=0;i<input.NumNo;i++){
		pic.No[i]=*(mvint+i);
	    }
	    pv=mxGetField(prhs[9],0,"pic");
	    mvint=mxGetPr(pv);
	    for(i=0;i<input.n;i++){
		pic.Y[i]=*(mvint+i)+input.bsize;
		pic.X[i]=*(mvint+i+input.n)+input.bsize;
		pic.I[i]=*(mvint+i+input.n+input.n);
	    }
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
	    mv.No[mv.P]=(double *)malloc(sizeof(double)*input.NumNo);
	    for(i=0;i<input.n;i++){
		mv.X[mv.P][i]=pic.X[i];
		mv.Y[mv.P][i]=pic.Y[i];
		mv.I[mv.P][i]=pic.I[i];
	    }
	    for(i=0;i<input.NumNo;i++){
		mv.No[mv.P][i]=pic.No[i];
	    }
	    mv.P++;
	}
    /*Run EM*/
	RunStep(&run,&pic,&pic0,PosFix,false,NULL,EveNoise);
    /*Judge break or not*/
        if(run.terminate=='Y')
	    break;
    }
    if(t==input.MaxIt)
        mexPrintf("\nReached maximal iteration number: %d iterations\n\n",t);
    else
	mexPrintf("\nFinished, converge after %d iterations\n\n",t);
    mv.T=t;
    WriteOutput(nlhs,plhs,&run,&pic,&mv,input.seed,input.bsize);
    FreeMem(&run,&input,&mv,&pic,&pic0);

    return;
}








