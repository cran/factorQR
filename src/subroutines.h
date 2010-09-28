#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include "vector.h"
//#include "rand.h"


void SWP( double **X, int k, int size);
void dinv(double **X, int size, double **X_inv);
void dcholdc(double **X, int size, double **L);
double ddet(double **X, int size, int give_log);
double checkFcn(double x, double pQuant);
/*void rinvGauss(double *normArray,int n, double mu,double lambda);*/
double rinvGauss2(double mu, double lambda);
void rgauss(double *normArray, int n, double mean, double sd);
double mhProb(double **invPhi, double **Lambda, double lambdaInt,
			double *OmegaI, double *invTildePsi, 
			double xBeta, int nFact, int nDimX, double *x, double y);
/*void rMVN(double *Sample, double *mean,double **Var,int size);
void rWish(double **Sample,double **S,int df,int size);
void dinv(double **X,int size, double **X_inv);
void SWP(double **X,int k,int size);      */

double TruncNorm(double lb, double ub, double mu, double var, int invcdf);
void rMVN(double *Sample, double *mean, double **inv_Var, int size);
void rWish(double **Sample, double **S, int df, int size);
