#ifndef EMfunctions_H
#define EMfunctions_H
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#define DoubleZero 1e-10
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define Zero 1e-1    /*Molecules below this intensity are vanished*/
#define MZero 1e-8

#define Expmin 32
#define Expinterval 1024 

struct Sparse{
    int n;
    short *X;
    short *Y;
    double *V;
};

struct RUN{
    int n;
    double sig;
    int bsize;
    int psfdecay;
    double Pzero;
    double Izero;
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

struct State{
    int n;
    double no;
    double *X;
    double *Y;
    double *I;
};

struct CPtable
{
    double C;
    double LogStepLambda;
    double LogStepMu;
    double LogMinLambda;
    double LogMaxLambda;
    double LogMinMu;
    double LogMaxMu;
    int NumMu;
    int NumLambda;
    double *MuX;
    double *LambdaX;
    double *LogMuX;
    double *LogLambdaX;
    double **CPtable;
};

void Sparse_malloc(struct Sparse *Sp, int n);
void Sparse_free(struct Sparse *Sp);
void State_malloc(struct State *St, int n);
void State_free(struct State *St);
void RUN_free(struct RUN *R);
void CPtable_malloc(struct CPtable *CP, int NumLambda, int NumMu);
void CPtalbe_free(struct CPtable *CP);
void CPtable_prepare(struct CPtable *CP);
double calibrated_poisson(double Lambda, double Mu, struct CPtable *CP);
void PSFGauss(struct Sparse *A, struct RUN *run, int L1, int U1, int L2, int U2, double x, double y, double sigma, int psfdecay);
void RunStep(struct RUN *run, struct State *pic, struct State *pic0);


#endif





