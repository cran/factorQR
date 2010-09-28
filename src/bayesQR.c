/*
 *  quantRegFact.c
 *  
 *
 *  Created by Lane Burgette on 2/23/10.
 *
 */

#include <string.h>
#include <stddef.h>
#include <stdio.h>      
#include <math.h>
#include <Rmath.h>
#include <R.h>
#include "vector.h"
#include "subroutines.h"
#include <R_ext/Lapack.h>
//#include "rand.h"



	
void cBayesQuantReg(int *piNObs,
					int *piNGen,
					int *piNDimX,
					double *pdX,
					double *pdY,
					double *pdXBeta, 
					double *piPQuant,
					int *piBurn,
					int *piThin,
					double *piC0,
					double *piD0,
					double *piB0q,
					double *piBeta,
					int *piPrint,
					double *pdStore){
	/*parameters from R */
	int nObs = *piNObs;
	int n_gen = *piNGen; 
	int nDimX = *piNDimX;
	double pQuant = *piPQuant;
	int thin = *piThin;
	double C0 = *piC0;
	double D0 = *piD0;
	
	double **X;
	double *Y;
	double *XBeta;
	
	double **B0q;
	double *vVec;
	double *uVec;	
	double *vtemp2;
	double **mtemp;
	double **betaVar;
	double *hatBeta;
	
	
	X = doubleMatrix(nObs, nDimX);
	Y = doubleArray(nObs);	
	XBeta = doubleArray(nObs);
	
	B0q = doubleMatrix(nDimX, nDimX);	
	vVec = doubleArray(nObs);
	uVec = doubleArray(nObs);
	vtemp2 = doubleArray(nDimX); //For sampling beta
	mtemp = doubleMatrix(nDimX, nDimX);
	betaVar = doubleMatrix(nDimX, nDimX);
	hatBeta = doubleArray(nDimX);
	
	
	int i, j,k, itemp, main_loop;
	int keep = 1, itempS = 0, itempP=ftrunc((double) n_gen/10), progress = 1;
	double tau, dtemp, dtemp2, dtemp3;
	
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		for(j=0; j<nObs; j++)
			X[j][i]=pdX[itemp++];
	
	double *beta;
	beta = doubleArray(nDimX);
		
		itemp = 0;
		for(i=0; i<nDimX; i++)
			beta[i] = piBeta[itemp++];
		
	itemp = 0;
	for(i=0; i<nObs; i++)
		XBeta[i] = pdXBeta[itemp++];
	
	
	itemp=0;
	for(i=0; i<nObs; i++)
		Y[i]=pdY[itemp++];
	
	

	
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		for(j=0; j<nDimX; j++)
			B0q[i][j] = piB0q[itemp++];
	
	
	
	
	GetRNGstate();
	
	for(main_loop=1; main_loop <= n_gen; main_loop++){
		
		if(main_loop != 1){
		for(i=0; i<nObs; i++){
			XBeta[i] = 0.0;
			for(j=0; j<nDimX; j++)
				XBeta[i] += X[i][j]*beta[j];
		}}
				
		
		
		
		dtemp=0;
		for(i=0; i<nObs; i++)
			dtemp += checkFcn(Y[i] - XBeta[i], pQuant);
		
		for(i=0; i<nDimX; i++)
			vtemp2[i] = 0;

		
		tau = (double) rgamma((C0 + (double) nObs), 1/(D0 + dtemp));
		
		
		dtemp2 = tau/(2* pQuant*(1-pQuant));
		for(i=0; i<nObs; i++){
			dtemp = 1/(pQuant * (1 - pQuant) * fabs(Y[i] - XBeta[i]));
			vVec[i] = rinvGauss2(dtemp, dtemp2);
			
		}
		
		
//		Rprintf("vVec: \n");
//		PdoubleArray(vVec, 10);		
		
		for(i=0; i<nObs; i++)
			uVec[i]=Y[i] - (1-2*pQuant)/(pQuant*(1-pQuant) * vVec[i]);
		
		
		for(i=0; i<nDimX; i++)
			for(j=0; j<nDimX; j++)
				mtemp[i][j] = 0;
		
		for(i=0; i<nDimX; i++)
			for(j=0; j<nDimX; j++)
				for(k=0; k<nObs; k++)
					mtemp[i][j] += X[k][i]*X[k][j] * vVec[k];
		
		
		for(i=0; i<nDimX; i++)
			for(j=0; j<nDimX; j++){
				mtemp[i][j] *= .5 * tau * pQuant * (1-pQuant); 
				mtemp[i][j] += B0q[i][j];
			}
		
		
		dinv(mtemp, nDimX, betaVar);
		
		
		
		for(i = 0; i<nDimX; i++){
			vtemp2[i] = 0; 
			for(j=0; j<nObs; j++)
				vtemp2[i] += .5*tau * pQuant *(1-pQuant) * X[j][i] * vVec[j] * uVec[j];
		}
		
		
		for(i=0; i<nDimX; i++)
			hatBeta[i] = 0;
		
		for(i=0; i<nDimX; i++)
			for(j=0; j<nDimX; j++)
				hatBeta[i] += betaVar[i][j] * vtemp2[j];
		
		

		rMVN(beta, hatBeta, betaVar, nDimX);
		
		for(j=0; j<nDimX; j++){
			if((beta[j] == beta[j]) == 0){
				error("NaN returned for beta");
			}}

		R_CheckUserInterrupt();
		if(main_loop > *piBurn) {
			if(keep == (thin+1)) {
				
				for(j = 0; j<nDimX; j++)
					pdStore[itempS++] = beta[j];
				
				pdStore[itempS++] = tau;
				
				keep = 1;
			}
			else {
				keep++;
			}
		}
		
		if(*piPrint) {
			if(main_loop == itempP) {
				Rprintf("%3d percent done.\n", progress * 10);
				itempP += ftrunc((double) n_gen/10); progress++;
				R_FlushConsole();
			}
		}
	}
	PutRNGstate();
	

	free(beta);
	FreeMatrix(X, nObs);
	free(Y);
	FreeMatrix(B0q, nDimX);
//	FreeMatrix(R0, nFact);
//	free(invPsi);
//	FreeMatrix(invPhi, nFact);
//	FreeMatrix(Lambda, nDimX+1);
//	FreeMatrix(Omega, nObs);
	free(XBeta);
	free(uVec);
	free(vVec);
//	free(invTildePsi);
//	FreeMatrix(preMult, nDimX);
//	FreeMatrix(OmegaCross, nFact);
//	free(indexFcn);
//	free(fixedRows);
//	free(fVec);
//	free(nMean);
//	free(vtemp);
	free(vtemp2);
	FreeMatrix(mtemp, nDimX);
	FreeMatrix(betaVar, nDimX);
	free(hatBeta);
	
	
	
}

