//#include <stdio.h>
#include <cstdio>

//extern "C" {
#include <R.h>
#include <Rmath.h>
//};

void truncNormImpute(double* Y, double* Yobserved, int* delta, double* mtotalfit, double* indiv_locations, double* sig, int nobs) 
{
	 // think carefully about the indexing
	 double xi, tmp, fitted_val, U;
	 for(int i=1; i <= nobs; i++) {
	     if(delta[i]==0) {
	         U = unif_rand();
	         fitted_val = mtotalfit[i] + indiv_locations[i];
	         xi = (Yobserved[i] - fitted_val)/sig[0];  
	         // xi is often way too large on the first few iterations !!!
	         // if the normalized difference is greater than 4, set xi=4.0
	         if(xi > 4) {
	            xi = 4.0;
	         }
	         tmp = U + (1.0 - U)*pnorm(xi, 0.0, 1.0, 1, 0); // + unif_rand()*pnorm(xi, 0.0,1.0, 0, 0);
             Y[i] = fitted_val + sig[0]*qnorm(tmp, 0.0, 1.0, 1, 0); 
         }
     }
}

void truncNormImpute_SP(double* Y, double* Yobserved, int* delta, double* mtotalfit, double sig, int nobs) 
{
	 // think carefully about the indexing
	 double xi, tmp, fitted_val, U;
	 for(int i=1; i <= nobs; i++) {
	     if(delta[i]==0) {
	         U = unif_rand();
	         fitted_val = mtotalfit[i];
	         xi = (Yobserved[i] - fitted_val)/sig;  
	         // xi is often way too large on the first few iterations !!!
	         // if the normalized difference is greater than 4, set xi=4.0
	         if(xi > 4) {
	            xi = 4.0;
	         }
	         tmp = U + (1.0 - U)*pnorm(xi, 0.0, 1.0, 1, 0); // + unif_rand()*pnorm(xi, 0.0,1.0, 0, 0);
             Y[i] = fitted_val + sig*qnorm(tmp, 0.0, 1.0, 1, 0); 
         }
     }
}

