#include"mex.h"
#include<stdio.h>
#include<stdlib.h>
/*Modified in 2015.06.19 by ZHANG Haowen
*/
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

void CalDist(struct MolList *M, struct Dist *D)
{
    int i,j,k;
    int n;
    n=M->n;
    for(i=0;i<n;i++){
	for(j=i;j<n;j++){
	    D->D[i][j]=(M->x[i]-M->x[j])*(M->x[i]-M->x[j])+(M->y[i]-M->y[j])*(M->y[i]-M->y[j]);
	    D->D[j][i]=D->D[i][j];
	}
    }
    return;
}

void DFS(struct Dist *D, double th, struct Group *G)
{
    int i,j,k;
    int n;
    n=D->n;
    char find[n];
    struct Node Link[n];
    struct Node *P,*TempP,*TempP2,*TempP3;
    for(i=0;i<n;i++){
	Link[i].next=NULL;
	for(j=0;j<n;j++){
	    if(i==j)
		continue;
	    if(D->D[i][j]<th*th){
		P=(struct Node *)malloc(sizeof(struct Node));
		P->value=j;
		TempP=Link[i].next;
		Link[i].next=P;
		P->next=TempP;
	    }
	}
    }
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
    struct Dist D;
    CreateDist(&D,n);
    CalDist(&ML,&D);
    struct Group G;
    DFS(&D,th,&G);
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
    DeleteDist(&D);
    return;
}









