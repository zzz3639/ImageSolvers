#include"mex.h"
#include"EMfunctions.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    struct CPtable CP;
    double *t;
    int NumLambda,NumMu;
    int i,j;
    mxArray *pv1,*pv2;
    pv1=mxGetField(prhs[nrhs-1],0,"C");
    t=mxGetPr(pv1);
    CP.C=(*t);
    pv1=mxGetField(prhs[nrhs-1],0,"MuX");
    i=mxGetM(pv1);
    j=mxGetN(pv1);
    NumMu=i*j;
    pv2=mxGetField(prhs[nrhs-1],0,"LambdaX");
    i=mxGetM(pv2);
    j=mxGetN(pv2);
    NumLambda=i*j;
    CPtable_malloc(&CP,NumLambda,NumMu);
    t=mxGetPr(pv1);
    for(i=0;i<NumMu;i++)
        CP.MuX[i]=*(t+i);
    t=mxGetPr(pv2);
    for(i=0;i<NumLambda;i++)
        CP.LambdaX[i]=*(t+i);
    
    pv1=mxGetField(prhs[nrhs-1],0,"CPtable");
    t=mxGetPr(pv1);
    for(i=0;i<NumLambda;i++)
        for(j=0;j<NumMu;j++)
	    CP.cptable[i][j]=*(t+j*NumLambda+i);

    double x,y;
    x=*(mxGetPr(prhs[0]));
    y=*(mxGetPr(prhs[1]));
    CPtable_prepare(&CP);

    mexPrintf("\n%f\n",calibrated_poisson(x,y,&CP));



    return;
}
