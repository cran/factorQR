
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Rmath.h>
#include <R.h> 
#include <R_ext/Lapack.h>
#include "vector.h"
//#include "rand.h"

/* Some of the following are taken from Imai and van Dyk's MNP package.
 See imai.princeton.edu/research/mnp.html 
 In particular, SWP, dinv, dcholdc, TruncNorm, rMVN, and rWish
 */

/*  The Sweep operator */
void SWP(
	 double **X,             /* The Matrix to work on */
	 int k,                  /* The row to sweep */
	 int size)               /* The dim. of X */
{
  int i,j;

  if (X[k][k] < 10e-20) 
    error("SWP: singular matrix.\n");
  else
    X[k][k]=-1/X[k][k];
  for(i=0;i<size;i++)
    if(i!=k){
      X[i][k]=-X[i][k]*X[k][k];
      X[k][i]=X[i][k];
    }
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      if(i!=k && j!=k)
	X[i][j]=X[i][j]+X[i][k]*X[k][j]/X[k][k];
  
}

double mhProb(
			double **invPhi,
			double **Lambda,
			double lambdaInt,
			double *OmegaI,
			double *invTildePsi,
			double xBeta,
			int nFact,
			int nDimX,
			double *x,
			double y)
{
	int i, j;
	double *vTemp = doubleArray(nFact);
	double dTemp1, dTemp2;
	for(j=0; j<nFact; j++){
		dTemp1=0.0;
		for(i = 0; i<nFact; i++)
			dTemp1 += OmegaI[i] * invPhi[j][i];
		vTemp[j] = dTemp1;
	}
	dTemp1 = 0;
	dTemp2 = 0;
	for(i=0; i<nFact; i++)
		dTemp1 += OmegaI[i] * vTemp[i]; // dTemp1 is omega_i' invPhi omega_i
	
	for(i=0; i<nFact; i++)
		dTemp2	+= Lambda[nDimX][i] * OmegaI[i];
	dTemp2 += OmegaI[0] * OmegaI[1] * lambdaInt;
	dTemp1 += pow((y - dTemp2 - xBeta), 2.0)*invTildePsi[nDimX]; // add (y-xBeta-fcn(Lambda))
	
	for(j=0; j<nDimX; j++){
		dTemp2 = 0;
		for(i=0; i<nFact; i++){
			dTemp2 += Lambda[j][i] * OmegaI[i];
		}
		dTemp2 -= x[j];
		dTemp2 *= dTemp2;
		dTemp2 *= invTildePsi[j];
		dTemp1 += dTemp2;
	}
	
	free(vTemp);
	return(-.5 * dTemp1);
	
	
	
}


/* inverting a matrix */
void dinv(double **X,
	  int	size,
	  double **X_inv)
{
  int i,j, k, errorM;
  double *pdInv = doubleArray(size*size);

  for (i = 0, j = 0; j < size; j++) 
    for (k = 0; k <= j; k++) 
      pdInv[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdInv, &errorM);
  if (!errorM) {
    F77_CALL(dpptri)("U", &size, pdInv, &errorM);
    if (errorM) {
      Rprintf("LAPACK dpptri failed, %d\n", errorM);
      error("Exiting from dinv().\n");
    }
  }
  else {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dinv().\n");
  }
  for (i = 0, j = 0; j < size; j++) {
    for (k = 0; k <= j; k++) {
      X_inv[j][k] = pdInv[i];
      X_inv[k][j] = pdInv[i++];
    }
  }

  free(pdInv);
}


/* Cholesky decomposition */
/* returns lower triangular matrix */
void dcholdc(double **X, int size, double **L)
{
  int i, j, k, errorM;
  double *pdTemp = doubleArray(size*size);

  for (j = 0, i = 0; j < size; j++) 
    for (k = 0; k <= j; k++) 
      pdTemp[i++] = X[k][j];
  F77_CALL(dpptrf)("U", &size, pdTemp, &errorM);
  if (errorM) {
    Rprintf("LAPACK dpptrf failed, %d\n", errorM);
    error("Exiting from dcholdc().\n");
  }
  for (j = 0, i = 0; j < size; j++) {
    for (k = 0; k < size; k++) {
      if(j<k)
	L[j][k] = 0.0;
      else
	L[j][k] = pdTemp[i++];
    }
  }

  free(pdTemp);
} 

/* calculate the determinant of the positive definite symmetric matrix
   using the Cholesky decomposition  */
double ddet(double **X, int size, int give_log)
{
  int i;
  double logdet=0.0;
  double **pdTemp = doubleMatrix(size, size);
  
  dcholdc(X, size, pdTemp);
  for(i = 0; i < size; i++)
    logdet += log(pdTemp[i][i]);

  FreeMatrix(pdTemp, size);
  if(give_log)
    return(2.0*logdet);
  else
    return(exp(2.0*logdet));
}

/* Check function */
double checkFcn(double x, double pQuant) {
	double ans = 0.0;
	if (x<0)
		ans =-1*x*(1-pQuant);
	if(x>0)
		ans = pQuant*x;
	return(ans);
}


double rinvGauss2(
			   double mu,
			   double lambda
			   )
{ 
	double b=0.5*mu/lambda;
	double a=mu*b;
	double c=4.0*mu*lambda;
	double d=mu*mu;
	
	if(mu == 0){
		return(0.0);
	}
	if((mu == mu) == 0){
		error("mu == NaN for inverse Gaussian.\n");
		return(0.0);
	}
	else{

	double chiSq = norm_rand();
			double u=unif_rand();
			double v=chiSq*chiSq;	// Chi-square with 1 df
			double x=mu+a*v-b*sqrt(c*v+d*v*v);	// Smallest root
			double resp=(u<(mu/(mu+x)))?x:d/x;	// Pick x with prob mu/(mu+x), else d/x;
	return(resp);
	}
}	

// random gaussian-- also from suppDists

void rgauss(
			double* normArray,
			int n,
			double mean,
			double sd
			)
{ 
	int i;
	
	for (i=0;i<n;i++)
		normArray[i]=rnorm(mean,sd);
	
	
}



/* Sample from a univariate truncated Normal distribution 
 (truncated both from above and below): choose either inverse cdf
 method or rejection sampling method. For rejection sampling, 
 if the range is too far from mu, it uses standard rejection
 sampling algorithm with exponential envelope function. */ 
double TruncNorm(
				 double lb,  /* lower bound */ 
				 double ub,  /* upper bound */
				 double mu,  /* mean */
				 double var, /* variance */
				 int invcdf  /* use inverse cdf method? */
				 ) {
	
	double z;
	double sigma = sqrt(var);
	double stlb = (lb-mu)/sigma;  /* standardized lower bound */
	double stub = (ub-mu)/sigma;  /* standardized upper bound */
	if(stlb >= stub)
		error("TruncNorm: lower bound is greater than upper bound\n");
	if (invcdf) {  /* inverse cdf method */
		z = qnorm(runif(pnorm(stlb, 0, 1, 1, 0), pnorm(stub, 0, 1, 1, 0)),
				  0, 1, 1, 0); 
	}
	else { /* rejection sampling method */
		double tol=2.0;
		double temp, M, u, exp_par;
		int flag=0;  /* 1 if stlb, stub <-tol */
		if(stub<=-tol){
			flag=1;
			temp=stub;
			stub=-stlb;
			stlb=-temp;
		}
		if(stlb>=tol){
			exp_par=stlb;
			while(pexp(stub,1/exp_par,1,0) - pexp(stlb,1/exp_par,1,0) < 0.000001) 
				exp_par/=2.0;
			if(dnorm(stlb,0,1,1) - dexp(stlb,1/exp_par,1) >=
			   dnorm(stub,0,1,1) - dexp(stub,1/exp_par,1)) 
				M=exp(dnorm(stlb,0,1,1) - dexp(stlb,1/exp_par,1));
			else
				M=exp(dnorm(stub,0,1,1) - dexp(stub,1/exp_par,1));
			do{ 
				u=unif_rand();
				z=-log(1-u*(pexp(stub,1/exp_par,1,0)-pexp(stlb,1/exp_par,1,0))
					   -pexp(stlb,1/exp_par,1,0))/exp_par;
			}while(unif_rand() > exp(dnorm(z,0,1,1)-dexp(z,1/exp_par,1))/M );  
			if(flag==1) z=-z;
		} 
		else{ 
			do z=norm_rand();
			while( z<stlb || z>stub ); 
		}
	}
	return(z*sigma + mu); 
}


/* Sample from the MVN dist */
void rMVN(                      
		  double *Sample,         /* Vector for the sample */
		  double *mean,           /* The vector of means */
		  double **Var,           /* The matrix Variance */
		  int size)               /* The dimension */
{
	int j,k;
	double **Model = doubleMatrix(size+1, size+1);
	double cond_mean;
    
	/* draw from mult. normal using SWP */
	for(j=1;j<=size;j++){       
		for(k=1;k<=size;k++)
			Model[j][k]=Var[j-1][k-1];
		Model[0][j]=mean[j-1];
		Model[j][0]=mean[j-1];
	}
	Model[0][0]=-1;
	Sample[0]=(double)norm_rand()*sqrt(Model[1][1])+Model[0][1];
	for(j=2;j<=size;j++){
		SWP(Model,j-1,size+1);
		cond_mean=Model[j][0];
		for(k=1;k<j;k++) cond_mean+=Sample[k-1]*Model[j][k];
		Sample[j-1]=(double)norm_rand()*sqrt(Model[j][j])+cond_mean;
	}
	
	FreeMatrix(Model,size+1);
}


/* Sample from a wish dist */
/* Odell, P. L. and Feiveson, A. H. ``A Numerical Procedure to Generate
 a Sample Covariance Matrix'' Journal of the American Statistical
 Association, Vol. 61, No. 313. (Mar., 1966), pp. 199-203. */

void rWish(                  
		   double **Sample,        /* The matrix with to hold the sample */
		   double **S,             /* The parameter */
		   int df,                 /* the degrees of freedom */
		   int size)               /* The dimension */
{
	int i,j,k;
	double *V = doubleArray(size);
	double **B = doubleMatrix(size, size);
	double **C = doubleMatrix(size, size);
	double **N = doubleMatrix(size, size);
	double **mtemp = doubleMatrix(size, size);
	
	for(i=0;i<size;i++) {
		V[i]=rchisq((double) df-i-1);
		B[i][i]=V[i];
		for(j=(i+1);j<size;j++)
			N[i][j]=norm_rand();
	}
	
	for(i=0;i<size;i++) {
		for(j=i;j<size;j++) {
			Sample[i][j]=0;
			Sample[j][i]=0;
			mtemp[i][j]=0;
			mtemp[j][i]=0;
			if(i==j) {
				if(i>0)
					for(k=0;k<j;k++)
						B[j][j]+=N[k][j]*N[k][j];
			}
			else { 
				B[i][j]=N[i][j]*sqrt(V[i]);
				if(i>0)
					for(k=0;k<i;k++)
						B[i][j]+=N[k][i]*N[k][j];
			}
			B[j][i]=B[i][j];
		}
	}
	
	dcholdc(S, size, C);
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			for(k=0;k<size;k++)
				mtemp[i][j]+=C[i][k]*B[k][j];
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			for(k=0;k<size;k++)
				Sample[i][j]+=mtemp[i][k]*C[j][k];
	
	free(V);
	FreeMatrix(B, size);
	FreeMatrix(C, size);
	FreeMatrix(N, size);
	FreeMatrix(mtemp, size);
}
