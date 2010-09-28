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


/*void tstRmvn(void){
	int N=50;
	double **cov;
	double *mean;
	cov = doubleMatrix(2,2);
	cov[0][0]=4;
	cov[0][1]=0;
	cov[1][1]=1;
	cov[1][0]=0;
	mean = doubleArray(2);
	mean[0]=1;
	mean[1]=2;
	double *hldVec;
	hldVec = doubleArray(2);
	GetRNGstate();
	for(int i=0; i<N; i++){
		rMVN(hldVec, mean, cov, 2);
		PdoubleArray(hldVec, 2);
	}
	PutRNGstate();
}*/

/*
void tstRWish(void){
	int N=50;
	double **cov;
	double *mean;
	int size=2;
	cov = doubleMatrix(2,2);
	cov[0][0]=4;
	cov[0][1]=0;
	cov[1][1]=1;
	cov[1][0]=0;
	double **hldMat;
	hldMat = doubleMatrix(2,2);
	GetRNGstate();
	for(int i=0; i<N; i++){
		rWish(hldMat, cov, 5, 2);
		PdoubleMatrix(hldMat, 2,2);
	}
	PutRNGstate();
}

*/


	
void cQuantRegFact(int *piNObs, // number of observations
					int *piNGen, // number of iterations to sample
					int *piNDimX, // dimension of r.h.s. of factor model
					int *piNFact, // number of factors
					int *piXBetLen, // length beta (not in factor)
					double *pdX,  // X matrix (related to factors)
					double *pdY,  //response
					double *pdNonFactorX, // X matrix not related to factors
					double *piPQuant,  // response quantile
					int *piNDichot,  // number of dichotomous variables
					int *piWhichDichot, //ones and zeros indicating dichot. variables
				    int *piDichotFcn, //places of dichotomous variables
					int *piIndexFcn,  // keeps track of active dimensions of Lambda
					int *piFixedRows, // indicates rows of elements fixed at 1 in Lambda
					int *piBurn,  // burn-in iterations
					int *piThin,  // thinning interval
					double *piC0, // start hyperparameters
					double *piD0,
					double *piBetaPsiZero,
					double *piAlphaPsiZero,
					double *piSig0,
					double *piMu0,
					double *piR0,
					double *piNu0,
					double *piB0q, // end hyperparameters
					double *piBeta, // start initial values
					double *piPhiZero,
					double *piInvPsiZero,
					double *piLambdaZero,
					double *piOmegaZero, // end initial values
					int *piPrint, // verbose 
				    int *piStoreOmega,  //store Omega values?
				    int *piInteract, // interact first two latent factors?
				    double *piInteractManifest, // for interaction w/ manifest indicators
				    int *piNManiInt, //number of manifest interactions
				    int *piWhichManiInd, //which components of omega_i have manifest interactions?
					double *pdStore){
	/*parameters from R */
	int nObs = *piNObs;
	int n_gen = *piNGen; 
	int nDimX = *piNDimX;
	int nFact = *piNFact;
	int nDichot = *piNDichot;
	int nDimBeta = *piXBetLen;
	int nFactBeta = nFact + nDimBeta;
	double pQuant = *piPQuant;
	int nManInt = *piNManiInt;
	int burn = *piBurn;
	int thin = *piThin;
	double C0 = *piC0;
	double D0 = *piD0;
	double betaPsiZero = *piBetaPsiZero;
	double alphaPsiZero = *piAlphaPsiZero;
	double sig0 = *piSig0;
	double mu0 = *piMu0;
	double nu0 = *piNu0;
	int nPrint = *piPrint;
	int interact = *piInteract;
	int nFactBetaInt = nFactBeta;
	if(interact == 1) nFactBetaInt++;
	nFactBetaInt += nManInt;
	
	
	
	double **X;
	double *Y;
	
	double **B0q;
	double **R0;
	double **Phi;
	double **invPhi;
	double *invPsi;
	double **Lambda;
	double **Omega;
	double *lambdaQOmega;	
	double *vVec;
	double *uVec;
	double *invTildePsi;
	double *xBeta;
	double **preMult;	
	double **OmegaCross;
	double **LambdaTilde;
	int *indexFcn;
	int *fixedRows;
	int *whichDichot;
	int *dichotFcn;
	double *fVec;
	double *nMean;
	double *vtemp;
	double *vtemp2;
	double **mtemp;
	double **lamQVar;
	double *hatLambdaQ;
	double *LambdaQ;
	double *LambdaQInt;
	double **interactManifest;
	int *whichManiInd;
	
	
	X = doubleMatrix(nObs, nDimX);
	double *Xi; //ith row of X
	Xi = doubleArray(nDimX);
	double *accHld;
	accHld = doubleArray(2);
	Y = doubleArray(nObs);	
	xBeta = doubleArray(nObs);
	B0q = doubleMatrix(nFactBetaInt, nFactBetaInt);
	R0 = doubleMatrix(nFact,nFact); /*Prior scale for inv Wish */
	Phi = doubleMatrix(nFact,nFact);
	invPhi = doubleMatrix(nFact, nFact);
	invPsi = doubleArray(nDimX);
	Lambda = doubleMatrix(nDimX+1,nFact);
	LambdaTilde = doubleMatrix(nDimX+1, nFact);
	whichDichot = intArray(nDimX);
	Omega = doubleMatrix(nObs, nFactBetaInt);
	if(nManInt)
	interactManifest = doubleMatrix(nObs, nManInt);
	else interactManifest = doubleMatrix(nObs, 1);
	lambdaQOmega = doubleArray(nObs);	
	vVec = doubleArray(nObs);
	uVec = doubleArray(nObs);
	invTildePsi = doubleArray(nDimX + 1); 
	preMult = doubleMatrix(nFact,nFact); //For sampling Omega rows
	OmegaCross = doubleMatrix(nFact, nFact);
	indexFcn = intArray(nDimX);
	fixedRows = intArray(nDimX);
	fVec = doubleArray(nDimX + 1);
	nMean = doubleArray(nFact); //holds mean for sampling Omega
	vtemp = doubleArray(nFact); //receives Omega sample
	vtemp2 = doubleArray(nFactBetaInt); //For sampling Lambda_q and beta
	mtemp = doubleMatrix(nFactBetaInt, nFactBetaInt);
	lamQVar = doubleMatrix(nFactBetaInt, nFactBetaInt);
	hatLambdaQ = doubleArray(nFactBetaInt);
	LambdaQ = doubleArray(nFactBetaInt);
	LambdaQInt = doubleArray(nManInt);

	
	
	int i, j,k,l, itemp, colHld, main_loop;
	int keep = 1, itempS = 0, itempP=ftrunc((double) n_gen/10), progress = 1;
	double tau, dtemp, dtemp2, dtemp3, lambdaInteract=0.0, accRat, yHld, xBetHld;
	
	//Load values from R:  

	itemp=0;
	if(nManInt > 0){
		whichManiInd = intArray(nManInt);
		for(i=0; i<nManInt; i++)
			whichManiInd[i] = piWhichManiInd[itemp++];
	}
	else {
		whichManiInd = intArray(1);
		whichManiInd[0] = 0;
	}
	
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		for(j=0; j<nObs; j++)
			X[j][i]=pdX[itemp++];
	
	double *beta;
	
	if(nDimBeta > 0){
		beta = doubleArray(nDimBeta);
		
		itemp = 0;
		for(i=0; i<nDimBeta; i++)
			beta[i] = piBeta[itemp++];
		
	}
	else{
		beta = doubleArray(2);
	}
	
	
	
	itemp=0;
	for(i=0; i<nObs; i++)
		Y[i]=pdY[itemp++];
	
	
	
	itemp=0;
	for(i=0; i<(nFact); i++)
		for(j=0; j<(nFact); j++)
			R0[i][j] = piR0[itemp++];
	

	itemp=0;
	for(i=0; i<nFactBetaInt; i++)
		for(j=0; j<nFactBetaInt; j++)
			B0q[i][j] = piB0q[itemp++];

	
	
	itemp=0;
	for(i=0; i<nFact; i++)
		for(j=0; j<nFact; j++)
			invPhi[i][j] = piPhiZero[itemp++];
	
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		invPsi[i] = piInvPsiZero[itemp++];
	
	
	itemp=0;
	for(j=0; j<nFact; j++)
		for(i=0; i<(nDimX + 1); i++)
			Lambda[i][j] = piLambdaZero[itemp++];
	
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		indexFcn[i] = piIndexFcn[itemp++];
	
	for(i=0; i<nDimX; i++)
		indexFcn[i] -= 1; //Change to C indexing from R
	
	itemp=0;
	for(i=0; i<nDimX; i++)
		whichDichot[i] = piWhichDichot[itemp++];
	

	
	itemp=0;
	for(i=0; i<nDimX; i++)
		fixedRows[i] = piFixedRows[itemp++];
	
	
	itemp=0;
	for(i=0; i<nObs; i++)
		for(j=0; j<nFact; j++)
			Omega[i][j] = piOmegaZero[itemp++];
	
	itemp=0;
	for(i=0; i<nObs; i++)
		for(j=0; j< nManInt; j++)
			interactManifest[i][j] = piInteractManifest[itemp++];
	
	if(nDimBeta > 0){
		itemp=0;
		for(i=0; i<nDimBeta; i++)
			for(j=0; j<nObs; j++)
				Omega[j][i+nFact]=pdNonFactorX[itemp++];	
	}
	
	if(interact == 1){
		for(i=0; i<nObs; i++)
			Omega[i][nFactBeta] = Omega[i][0] * Omega[i][1];
	}
	
	for(i=0; i<nObs; i++)
		xBeta[i] = 0.0;
	
	

		double **xForDichot;
		double *muXDichot;
//vector of integers w/ places of dichotomous x's
//	int *dichotInd; //vector of 1's and 0's for dichotomous x's 
	
	
	GetRNGstate();	

	if(nDichot > 0){	
		dichotFcn = intArray(nDichot);
		xForDichot= doubleMatrix(nObs, nDichot);
		muXDichot = doubleArray(nDichot);
//		dichotInd = intArray(nDichot);
		itemp=0;
		for(i=0; i<nDichot; i++)
			dichotFcn[itemp++] = piDichotFcn[i] - 1;
		
		for(i=0; i<nDichot; i++)
			muXDichot[i] = 0.0;
		
		for(i=0; i<nObs; i++)
			for(j=0; j<nDichot; j++)
				xForDichot[i][j] = X[i][dichotFcn[j]];
		
		for(i=0; i<nObs; i++)
			for(j=0; j< nDichot; j++){
				if(xForDichot[i][j] > 0)
					X[i][dichotFcn[j]] = TruncNorm(0.0, 1000.0, 0.0, 1.0, 0);
				else
					X[i][dichotFcn[j]] = TruncNorm(-1000.0, 0.0,0.0, 1.0, 0);
				}

		
	}
	
	else{
		xForDichot=doubleMatrix(2,2);
		muXDichot=doubleArray(2);
		dichotFcn = intArray(2);
		
	}
	

	

	
	for(main_loop=1; main_loop <= n_gen; main_loop++){
		
		for(i=0; i<nObs; i++){
			lambdaQOmega[i] = 0.0;
			for(j=0; j<nFact; j++)
				lambdaQOmega[i] += Omega[i][j]*Lambda[nDimX][j];
		}
		
		
// last column of Omega is interaction if needed:		
		if(interact == 1){
			for(i = 0; i<nObs; i++)
				lambdaQOmega[i] += Omega[i][0] * Omega[i][1] * lambdaInteract;
		}
		
//		Rprintf("lambdaQOmega: \n");
//		PdoubleArray(lambdaQOmega, 10);
		
		if(nDimBeta > 0){
			for(i=0; i<nObs; i++)
				for(j=0; j<nDimBeta; j++)
					lambdaQOmega[i] += Omega[i][j+nFact] * beta[j];
		}
		
		

		if(nManInt > 0){
			for(i=0; i<nObs; i++)
				for(k=0; k<nManInt; k++)
					lambdaQOmega[i] += LambdaQInt[j] * Omega[i][whichManiInd[k]]*interactManifest[i][k];
		}
		
		dtemp=0;
		for(i=0; i<nObs; i++)
			dtemp += checkFcn(Y[i] - lambdaQOmega[i], pQuant);
		
		
		
// draw inverse scale of response:		

		tau = (double) rgamma((C0 + (double) nObs), 1/(D0 + dtemp));
		
		
		dtemp2 = tau/(2* pQuant*(1-pQuant));
		for(i=0; i<nObs; i++){
			dtemp = 1/(pQuant * (1 - pQuant) * fabs(Y[i] - lambdaQOmega[i]));
			vVec[i] = rinvGauss2(dtemp, dtemp2);
			
		}
		
		
//		Rprintf("vVec: \n");
//		PdoubleArray(vVec, 10);		
		
		for(i=0; i<nObs; i++)
			uVec[i]=Y[i] - (1-2*pQuant)/(pQuant*(1-pQuant) * vVec[i]);
		
//		Rprintf("uVec: \n");
//		PdoubleArray(uVec, 10);		
		
		
		/* Sample Omega */
		for(l=0; l<nObs; l++){	
			if(nManInt > 0){
				for(j=0; j<nFact; j++)
					for(k=0; k<(nDimX + 1); k++)
						LambdaTilde[k][j] = Lambda[k][j];
				
				for(i=0; i<nManInt; i++)
					LambdaTilde[nDimX][whichManiInd[i]] += interactManifest[l][i] * LambdaQInt[i];
			}
			else {
				if (l == 0) {
					for(j=0; j<nFact; j++)
						for(k=0; k<(nDimX + 1); k++)
							LambdaTilde[k][j] = Lambda[k][j];
				}
		
			}

			
			for(i=0; i<nFact; i++)
				for(j=0; j<(nFact); j++){
					preMult[i][j]=0;
					for(k=0; k<(nDimX + 1); k++){
						invTildePsi[k] = (k < nDimX ? invPsi[k] : .5*tau*pQuant*(1-pQuant)*vVec[l]); 
						preMult[i][j] += LambdaTilde[k][i] * invTildePsi[k] * LambdaTilde[k][j];
					}
				}
			

			
			for(i=0; i<nFact; i++)
				for(j=0; j<nFact; j++)
					preMult[i][j] += invPhi[i][j];
			
			/*			void dinv(double **X, int size, double **X_inv);*/
			
//			Rprintf("preMult \n");
//			PdoubleMatrix(preMult, nFact, nFact);
			
			dinv(preMult, nFact, mtemp);		
			
//			Rprintf("mtemp \n");
//			PdoubleMatrix(mtemp, nFact, nFact);

			
			
			for(i=0; i< nFact; i++)
				fVec[i] = 0.0;
			
			for(j=0; j<nFact; j++)
				for(k=0; k<(nDimX + 1); k++){
					/*fVec[j] += Lambda[k][j] * invTildePsi[k] * (k<nDimX ? X[l][k] : uVec[k]); */
					fVec[j] += LambdaTilde[k][j] * invTildePsi[k] * (k<nDimX ? X[l][k] : uVec[l]);
				}
			
//			Rprintf("fVec \n");
//			PdoubleArray(fVec, nFact);
			
			
			for(i=0; i<nFact; i++){
				nMean[i] = 0;
				for(j=0; j<nFact; j++)
					nMean[i] += mtemp[i][j]*fVec[j];
			}
			
			
//			Rprintf("nMean \n");
//			PdoubleArray(nMean, nFact);
			 
			 /*			void rMVN(double *Sample, double *mean, double **inv_Var, int size);
			 It should be (Sample, mean, Var, size), not inv_Var*/
			/*rMVN(vtemp, nMean, preMult, nFact);*/
			if(interact == 1){
				if(nDimBeta > 0){
					for(i=0; i<nObs; i++){
						xBeta[i] = 0;
						for(j=0; j<nDimBeta; j++)
							xBeta[i] += Omega[i][j+nFact] * beta[j];
					}
				}
				
				for(i=0; i<nDimX; i++)
					Xi[i] = X[l][i];
				for(i=0; i<nFact; i++)
					nMean[i] = Omega[l][i];
				rMVN(vtemp, nMean, mtemp, nFact);
//								Rprintf("before accRat \n");
//				R_FlushConsole();
//				yHld = uVec[l];
//				xBetHld = xBeta[l];
				accRat = mhProb(invPhi, Lambda, lambdaInteract, vtemp, invTildePsi, xBeta[l], 
					   nFact, nDimX, Xi, uVec[l]);
				accRat = accRat - mhProb(invPhi, Lambda, lambdaInteract, nMean, invTildePsi, xBeta[l], 
								nFact, nDimX, Xi, uVec[l]);				
				//accRat = -2;
				//Rprintf("after accRat; log accRat %f \n", accRat);
				accRat = exp(accRat);
				if(unif_rand() <= accRat)
					for(i=0; i<nFact; i++){
						Omega[l][i] = vtemp[i];}
				
				Omega[l][nFactBeta] = Omega[l][0]*Omega[l][1];
			}
			else{
				rMVN(vtemp, nMean, mtemp, nFact);
//			PdoubleArray(vtemp, nFact);
				for(i=0; i<nFact; i++)
					Omega[l][i] = vtemp[i];
			}
			
		}
		
//		PdoubleMatrix(Omega, 10, nFactBetaInt);
		
// Sample Lambda_q (and beta if needed)
		
//		Rprintf("made it out of Omega \n");
//		R_FlushConsole();	
		
		if(nManInt > 0)
			for(i=0; i<nObs; i++)
				for(j=0; j<nManInt; j++)
					Omega[i][j+nFactBeta] = Omega[i][whichManiInd[j]]*interactManifest[i][j];
		
		for(i=0; i<nFactBetaInt; i++)
			for(j=0; j<nFactBetaInt; j++)
				mtemp[i][j] = 0;
		
		for(i=0; i<nFactBetaInt; i++)
			for(j=0; j<nFactBetaInt; j++)
				for(k=0; k<nObs; k++)
					mtemp[i][j] += Omega[k][i]*Omega[k][j] * vVec[k];
		
		
		for(i=0; i<nFactBetaInt; i++)
			for(j=0; j<nFactBetaInt; j++){
				mtemp[i][j] *= .5 * tau * pQuant * (1-pQuant); 
				mtemp[i][j] += B0q[i][j];
			}
		
//		Rprintf("mtemp \n");
//		PdoubleMatrix(mtemp, nFactBeta, nFactBeta);
		
		dinv(mtemp, nFactBetaInt, lamQVar);
		
//		PdoubleMatrix(mtemp, nFactBetaInt,nFactBetaInt);
//		R_FlushConsole();
		
//		Rprintf("lamQVar \n");
//		PdoubleMatrix(lamQVar, nFactBeta, nFactBeta);
		
		
		for(i = 0; i<nFactBetaInt; i++){
			vtemp2[i] = 0; 
			for(j=0; j<nObs; j++)
				vtemp2[i] += .5*tau * pQuant *(1-pQuant) * Omega[j][i] * vVec[j] * uVec[j];
		}
		
//		Rprintf("vtemp2 \n");
//		PdoubleArray(vtemp2, nFactBeta);
		
		for(i=0; i<nFactBetaInt; i++)
			hatLambdaQ[i] = 0;
		
		for(i=0; i<nFactBetaInt; i++)
			for(j=0; j<nFactBetaInt; j++)
				hatLambdaQ[i] += lamQVar[i][j] * vtemp2[j];
		
		

		rMVN(LambdaQ, hatLambdaQ, lamQVar, nFactBetaInt); 

		
		for(i=0; i<nFact; i++)
			Lambda[nDimX][i] = LambdaQ[i];
		
		if(nDimBeta > 0){
			for(i=0; i<nDimBeta; i++)
				beta[i] = LambdaQ[nFact + i];
		}
		
		if(interact == 1)
			lambdaInteract = LambdaQ[nFactBeta];
		
		if(nManInt > 0)
			for(i=0; i< nManInt; i++){
				LambdaQInt[i] = LambdaQ[nFactBeta + i];
			}
		
//		PdoubleArray(LambdaQ, nFactBetaInt);
		
		
		for(i=0; i<nDimX; i++){
			colHld = indexFcn[i];
			if(fixedRows[i] == 1){
				if(whichDichot[i] == 1){
					invPsi[i] = 1.0;
				}
				else{
				dtemp = 0;
				for(j=0; j<nObs; j++)
					dtemp += pow(X[j][i] - Omega[j][colHld], 2.0);
				invPsi[i] = (double) rgamma((double) nObs/2 + alphaPsiZero, 1/(betaPsiZero + .5*dtemp));
				}
			}
			else {
				dtemp = 0;
				dtemp2 = 0;
				for(j=0; j<nObs; j++){
					dtemp += Omega[j][colHld] * X[j][i];
					dtemp2 += pow(Omega[j][colHld], 2.0);
				}
				dtemp3 = (dtemp + sig0 * mu0)/(sig0 + dtemp2);
				dtemp = 0;
				dtemp2 = 0;
				for(j=0; j<nObs; j++){
					dtemp += pow(X[j][i] - dtemp3 * Omega[j][colHld],2.0);
					dtemp2 += pow(Omega[j][colHld],2.0);
				}
				dtemp += pow(dtemp3 - mu0, 2.0) *sig0;
				dtemp *= .5;
				dtemp += betaPsiZero;
				if(whichDichot[i] == 1){
					invPsi[i] = 1.0;
				}
				else{
				invPsi[i] = (double) rgamma((double)nObs/2.0 + alphaPsiZero, 1/dtemp);
				}
				Lambda[i][colHld] = dtemp3 + norm_rand() / sqrt(invPsi[i]*(sig0 + dtemp2));
			}
			
		}

		
		
		for(j=0; j<nFact; j++)
			for(k=0; k<nFact; k++)
				OmegaCross[j][k] = R0[j][k];
		
		for(i=0; i<nFact; i++)
			for(j=0; j<nFact; j++)
				for(k=0; k<nObs; k++)
					OmegaCross[i][j] += Omega[k][i] * Omega[k][j];
		
		dinv(OmegaCross, nFact, Phi);
		rWish(invPhi, Phi, nObs + nu0, nFact);
		dinv(invPhi, nFact, Phi);
		
		if(nDichot > 0){
			for(i=0; i<nObs; i++)
				for(j=0; j< nDichot; j++){
					dtemp = 0.0;
					for(k=0; k<nFact; k++){
						dtemp += Omega[i][k] * Lambda[dichotFcn[j]][k];
					}
					if(xForDichot[i][j] > 0)
						X[i][dichotFcn[j]] = TruncNorm(0.0, 1000.0, muXDichot[j] + dtemp, 1.0, 0);
					else
						X[i][dichotFcn[j]] = TruncNorm(-1000.0, 0.0,muXDichot[j] + dtemp, 1.0, 0);
				}	
			
			for(i=0; i<nDichot; i++){
				itemp = dichotFcn[i];
				dtemp = 0.0;
				for(j=0; j< nObs; j++)
					for(k = 0; k<nFact; k++)
						dtemp += Lambda[itemp][k] * Omega[j][k];
				dtemp2 = 1.0 + (double) nObs * invPsi[itemp];
				
//				PintArray(dichotFcn, 1);
//				Rprintf("%f dtemp2 \n", dtemp2);
//				Rprintf("%f dtemp \n\n", dtemp);
//				PdoubleMatrix(X, 20, nDimX);
//				PdoubleMatrix(Omega, 10, 1);
//				PdoubleMatrix(xForDichot, 20, 1);
				
						
//				muXDichot[i] = ((double) nObs * invPsi[itemp] * dtemp)/dtemp2 +  norm_rand()/sqrt(dtemp2);
				muXDichot[i] = (invPsi[itemp] * dtemp)/dtemp2 +  norm_rand()/sqrt(dtemp2);
				
			}
			
			//change to X-muXDichot since all draws except mu use that instead of X
			
			for(i=0; i<nObs; i++)
				for(j=0; j<nDichot; j++)
					X[i][dichotFcn[j]] -= muXDichot[j];
				
		}
			
			
			
		
		R_CheckUserInterrupt();
		if(main_loop > *piBurn) {
			if(keep == (thin+1)) {
				for(j=0; j<nFact; j++)
					pdStore[itempS++] = Lambda[nDimX][j];
				
				if(interact == 1)
					pdStore[itempS++] = lambdaInteract;				
				
				if(nManInt > 0)
					for(j=0; j<nManInt; j++)
						pdStore[itempS++] = LambdaQInt[j];				
				
				for(j=0; j<nDimX; j++)
					pdStore[itempS++] = Lambda[j][indexFcn[j]];
				
				if(nDimBeta > 0)
					for(j = 0; j<nDimBeta; j++)
						pdStore[itempS++] = beta[j];
				

				
				for(j=0; j<nFact; j++)
					for(k=0; k<nFact; k++)
						if(j<=k)
							pdStore[itempS++] = Phi[j][k];
				for(j=0; j<nDimX; j++)
					pdStore[itempS++] = 1.0/invPsi[j];
				

				
				pdStore[itempS++] = tau;
				
				if(nDichot > 0)
					for(j =0; j<nDichot; j++)
						pdStore[itempS++] = muXDichot[j];
				
				if(*piStoreOmega){
					for(j = 0; j< nObs; j++){
						for(k=0; k<nFact; k++){
							pdStore[itempS++] = Omega[j][k];
						}
					}
				}
				
				keep = 1;
			}
			else {
				keep++;
			}
		}
		
		if(nPrint) {
			if(main_loop == itempP) {
				Rprintf("%3d percent done.\n", progress * 10);
				itempP += ftrunc((double) n_gen/10); progress++;
				R_FlushConsole();
			}
		}
	}
	PutRNGstate();
	

	free(beta);
	free(muXDichot);
	free(dichotFcn);
	if(nDichot > 0)
		FreeMatrix(xForDichot, nObs);
	else 
		FreeMatrix(xForDichot, 2);
	FreeMatrix(X, nObs);
	free(Y);
	free(xBeta);
	FreeMatrix(B0q, nFactBetaInt);
	FreeMatrix(R0, nFact);
	free(invPsi);
	FreeMatrix(invPhi, nFact);
	FreeMatrix(Lambda, nDimX+1);
	FreeMatrix(LambdaTilde, nDimX+1);
	FreeMatrix(Omega, nObs);
	free(lambdaQOmega);
	free(uVec);
	free(vVec);
	free(invTildePsi);
	FreeMatrix(preMult, nFact);
	FreeMatrix(OmegaCross, nFact);
	free(indexFcn);
	free(fixedRows);
	free(fVec);
	free(nMean);
	free(vtemp);
	free(vtemp2);
	FreeMatrix(mtemp, nFactBetaInt);
	FreeMatrix(lamQVar, nFactBetaInt);
	FreeMatrix(interactManifest, nObs);
	free(Xi);
	free(LambdaQInt);
	
	
}

