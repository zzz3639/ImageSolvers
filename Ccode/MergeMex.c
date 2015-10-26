#include"mex.h"
#include<stdio.h>
#include<stdlib.h>
/*Modified in 2015.10.26 by ZHANG Haowen
*/
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
struct MolList
{
    int n;
    double *x;
    double *y;
};

struct Dist
{
    int n;
    double **D;
};

struct Node
{
    int value;
    struct Node *next;
};

struct Group
{
    struct Node *node;
    struct Group *next;
};

void CreateMolList(struct MolList *M, int n)
{
    M->n=n;
    M->x=(double *)malloc(sizeof(double)*n);
    M->y=(double *)malloc(sizeof(double)*n);
    return;
}

void DeleteMolList(struct MolList *M)
{
    free(M->x);
    free(M->y);
    return;
}

void CreateDist(struct Dist *D, int n)
{
    int i;
    D->n=n;
    D->D=(double **)malloc(sizeof(double *)*n);
    for(i=0;i<n;i++){
	D->D[i]=(double *)malloc(sizeof(double)*n);
    }
    return;
}

void DeleteDist(struct Dist *D)
{
    int i;
    for(i=0;i<D->n;i++){
	free(D->D[i]);
    }
    free(D->D);
    return;
}

struct Node* CalDist(struct MolList *M, double th)
{
    struct Node *Link;
    int i,j,k;
    int n;
    n=M->n;
    /*Calculate the minimal and maximal value of the Molecule list in both dimensions*/
    double Dzero=1e-10;
    double Xmin=1e10,Xmax=-1e10,Ymin=1e10,Ymax=-1e10;
    for(i=0;i<M->n;i++){
	if(M->x[i]<Xmin)
	    Xmin=M->x[i];
	if(M->x[i]>Xmax)
	    Xmax=M->x[i];
	if(M->y[i]<Ymin)
	    Ymin=M->y[i];
	if(M->y[i]>Ymax)
	    Ymax=M->y[i];
    }
    /*Indexing the molecules*/
    double Interval=MAX(th,0.3);
    int nx,ny;
    nx=(int)((Xmax-Xmin+Dzero)/Interval)+1;
    ny=(int)((Ymax-Ymin+Dzero)/Interval)+1;
    struct Node **PositionIndex;
    PositionIndex=(struct Node **)malloc(sizeof(struct Node *)*nx);
    for(i=0;i<nx;i++){
	PositionIndex[i]=(struct Node *)malloc(sizeof(struct Node)*ny);
    }
    for(i=0;i<nx;i++){
	for(j=0;j<ny;j++){
	    PositionIndex[i][j].next=NULL;
	}
    }
    int ix,iy;
    struct Node *p,*tempp;
    for(i=0;i<n;i++){
	ix=(int)((M->x[i]-Xmin+Dzero)/Interval);
	iy=(int)((M->y[i]-Ymin+Dzero)/Interval);
	p=(struct Node *)malloc(sizeof(struct Node));
	p->value=i;
	p->next=PositionIndex[ix][iy].next;
	PositionIndex[ix][iy].next=p;
    }

    /*Make the linking list*/
    int iu,ju,id,jd;
    double th2=th*th;
    Link=(struct Node *)malloc(sizeof(struct Node)*n);
    for(i=0;i<n;i++){
	Link[i].next=NULL;
	ix=(int)((M->x[i]-Xmin+Dzero)/Interval);
	iy=(int)((M->y[i]-Ymin+Dzero)/Interval);
	iu=MAX(0,ix-1);
	ju=MAX(0,iy-1);
	id=MIN(nx-1,ix+1);
	jd=MIN(ny-1,iy+1);
	for(k=iu;k<=id;k++){
	    for(j=ju;j<=jd;j++){
		tempp=PositionIndex[k][j].next;
		while(tempp!=NULL){
		    if(tempp->value==i){
			tempp=tempp->next;
			continue;
		    }
		    if((M->x[i]-M->x[tempp->value])*(M->x[i]-M->x[tempp->value])+(M->y[i]-M->y[tempp->value])*(M->y[i]-M->y[tempp->value])<th2){
			p=(struct Node *)malloc(sizeof(struct Node));
			p->value=tempp->value;
			p->next=Link[i].next;
			Link[i].next=p;
		    }
		    tempp=tempp->next;
		}
	    }
	}
    }

    /*free the indexing table*/
    for(i=0;i<nx;i++){
	for(j=0;j<ny;j++){
	    tempp=PositionIndex[i][j].next;
	    while(tempp!=NULL){
		p=tempp;
		tempp=tempp->next;
		free(p);
	    }
	}
	free(PositionIndex[i]);
    }
    free(PositionIndex);
    
    return Link;
}

void DFS(struct Node *Link, int n, struct Group *G)
{
    int i,j,k;
    char find[n];
    struct Node *P,*TempP,*TempP2,*TempP3;
    for(i=0;i<n;i++){
	find[i]=0;
    }
    struct Node stack;
    G->next=NULL;
    struct Group *PG,*TempPG;
    while(1){
	for(k=0;k<n;k++){
	    if(find[k]==0)
		break;
	}
	if(k==n)
	    break;
	PG=(struct Group *)malloc(sizeof(struct Group));
	PG->node=NULL;
	TempPG=G->next;
	G->next=PG;
	PG->next=TempPG;
	/*eliminate stack;*/
	P=(struct Node *)malloc(sizeof(struct Node));
        stack.next=P;
	P->value=k;
	P->next=NULL;
        find[k]=1;
	while(stack.next!=NULL){
	    P=stack.next;
	    stack.next=stack.next->next;
            
            
	    TempP=Link[P->value].next;
	    while(TempP!=NULL){
		if(find[TempP->value]==1){
		    TempP=TempP->next;
		    continue;
		}
                find[TempP->value]=1;
		TempP2=(struct Node *)malloc(sizeof(struct Node));
                TempP2->value=TempP->value;
		TempP3=stack.next;
		stack.next=TempP2;
		TempP2->next=TempP3;
		TempP=TempP->next;
	    }

	    TempP=PG->node;
	    PG->node=P;
	    P->next=TempP;
	}
    }
    return;

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /*Reading the data*/
    if(nrhs!=2){
        mexPrintf("\nM=MergeMex(pic,th)\n");
	return;
    }
    double th;
    double *pointread;
    pointread=mxGetPr(prhs[1]);
    th=*pointread;
    int n;
    n=mxGetM(prhs[0]);
    pointread=mxGetPr(prhs[0]);
    struct MolList ML;
    CreateMolList(&ML,n);
    int i,j,k;
    for(i=0;i<n;i++){
	*(ML.x+i)=*(pointread+i);
	*(ML.y+i)=*(pointread+n+i);
    }
    /*Group the molecules with DFS*/
    struct Node *Link;
    Link=CalDist(&ML,th);
    struct Group G;
    DFS(Link,ML.n,&G);
    /*Grouping data processing and output*/
    struct Group *TempPG;
    struct Node *TempP;
    int ng=0;
    TempPG=G.next;
    while(TempPG!=NULL){
        ng+=1;
	TempPG=TempPG->next;
    }
    double *pointwrite;
    plhs[0]=mxCreateNumericMatrix(n,4,mxDOUBLE_CLASS,mxREAL);
    pointwrite=mxGetPr(plhs[0]);
    TempPG=G.next;
    i=0;
    k=0;
    while(TempPG!=NULL){
	TempP=TempPG->node;
	while(TempP!=NULL){
	    *(pointwrite+i)=ML.x[TempP->value];
	    *(pointwrite+i+n)=ML.y[TempP->value];
	    *(pointwrite+i+2*n)=TempP->value;
	    *(pointwrite+i+3*n)=k;
	    i+=1;
	    TempP=TempP->next;
	}
	k+=1;
	TempPG=TempPG->next;
    }
    DeleteMolList(&ML);
    /*free linking list*/
    struct Node *TempP2;
    for(i=0;i<n;i++){
	TempP2=Link[i].next;
	while(TempP2!=NULL){
            TempP=TempP2;
	    TempP2=TempP2->next;
	    free(TempP);
	}
    }
    free(Link);
    /*free grouping memory*/
    TempPG=G.next;
    while(TempPG!=NULL){
	TempP2=TempPG->node;
	while(TempP2!=NULL){
            TempP=TempP2;
	    TempP2=TempP2->next;
	    free(TempP);
	}
	G.next=TempPG;
	TempPG=TempPG->next;
	free(G.next);
    }
    return;
}









