#ifndef GUARD_Mixdev
#define GUARD_Mixdev

#include "Lib.h"
#include <numeric>

class Mixdev
{
public:

	// initialization
    Mixdev(): psi1(2.0), psi2(0.1) {}
	~Mixdev() {}
	// public methods
	void updateLabels(int* labels, double* mix_prop, double* locations, double* rr, double* sig, int nobs, int nclust) {
	    // Also need residuals, locations and sig!!!
	    // function to update cluster labels
	    // Remember that rr has nobs + 1 elements and starts from index 1.
	    double pp, ptot, pcum, uu;
	    int tmp;
	    int count;
	    for(int i=0; i < nobs; i++) {
	        ptot = 0.0;
	        count = 1;
	        for(int h=0; h < nclust; h++) {
	             ptot += mix_prop[h]*dnorm(rr[i+1], locations[h], sig[0], 0);
	        }
	        // sample uniformly from [0,ptot]
	        pcum = 0.0;
	        uu = ptot*unif_rand();
	        for(int h=0; h < nclust; h++) {
	             pcum += mix_prop[h]*dnorm(rr[i+1], locations[h], sig[0], 0);
	             if(uu < pcum) {
	                break;
	             }
	             count = count + 1;
	        }
	        labels[i] = count;
	     }   
	}
	void tabCounts(int*nn, int* labels, int nobs, int nclust) {
	     for(int h=0; h < nclust; h++) {
             nn[h] = std::count(labels, labels+nobs, h+1);
            // Rprintf("Cluster %d: %d \n", h+1, nn[h]);
         }
      //   Rprintf("Cluster 1: %d \n", nn[0]);
      //   Rprintf("Cluster 2: %d \n", nn[1]);
 	}
    void updateMix(double* mix_prop, double* mass, int* labels, int* nn, int nobs, int nclust)
	{    
	     // function to update mixture proportions and mass parameter.
	     double Vold, Vnew, log_mix_prop, Vcum;
	     double shape1, shape2;
	     double gam_shape, gam_scale, tmp;
	     int ncum;
         
         ncum = std::accumulate(nn + 1,nn + nclust,0);
         shape1 = nn[0] + 1.0;
         shape2 = mass[0] + ncum;
         Vold = rbeta(shape1, shape2);
         log_mix_prop = log(Vold);
         Vcum = log(1 - Vold);
         mix_prop[0] = Vold;
         for(int h=1; h < nclust - 1; h++) {
             ncum = std::accumulate(nn + h + 1,nn + nclust,0);
             shape1 = nn[h] + 1.0;
             shape2 = mass[0] + ncum;
             Vnew = rbeta(shape1, shape2);
             tmp = (Vnew/Vold)*(1 - Vold);
             log_mix_prop = log(tmp) + log_mix_prop;

             //log_mix_prop = log(Vnew) - log(Vold) + log(1.0 - Vold) + log_mix_prop;
             Vcum += log1p(-Vnew);
             Vold = Vnew;
             mix_prop[h] = exp(log_mix_prop);
         }
         mix_prop[nclust - 1] = 1.0 - std::accumulate(mix_prop,mix_prop + nclust - 1,0.0);
         
         gam_shape = psi1 + nclust - 1;
         gam_scale = 1/(psi2 - Vcum);
         mass[0] = rgamma(gam_shape, gam_scale);
         //Rprintf("Gamma draw: %f \n", Vcum);
 	}
    void updateLocations(double* locations, double* mix_prop, double* rr, int* labels, int* nn, double* sig, double prior_sigsq, int nobs, int nclust)
	{
	     double post_mean, post_var, post_sd, sigsq, dummy_bool, clust_sum, wts, muG;
	     sigsq = sig[0]*sig[0];
	     for(int h=0; h < nclust; h++) {
	         clust_sum = 0.0;
	         for(int i=0; i < nobs; i++) {
	             if(labels[i] == h+1) {
	                 dummy_bool = 1.0;
	             }
	             else {
	                 dummy_bool = 0.0;
	             }
	             clust_sum += dummy_bool*rr[i+1]; 
	         }
	         wts = prior_sigsq/(prior_sigsq*nn[h] + sigsq);  
	         post_mean = wts*clust_sum;
             post_var = sigsq*wts;
             post_sd = sqrt(post_var);
             locations[h] = rnorm(post_mean, post_sd);
         }    
         // Normalize the locations
         muG = 0.0;
         for(int h=0; h < nclust; h++) {
             muG = muG + mix_prop[h]*locations[h];
         }      
         for(int h=0; h < nclust; h++) {
             locations[h] = locations[h] - muG;
         }
	}
	void getIndivLocations(double* indiv_locations, double* locations, int* labels, int nobs, int nclust)
	{
	      double dum_location;
	      // index of indiv_locations starts from 1.
	      for(int i=0; i < nobs; i++) {
	           // find cluster observation i belongs to
	           for(int h=0; h < nclust; h++) {
	               if(labels[i] == h+1) {
	                    dum_location = locations[h];
	                    break;
	               }
	           }
	           indiv_locations[i+1] = dum_location;
	       }
	}
	void updateSigma(int* labels, double* locations, double* rr, double* sig, double kappa, int sigdf, int nobs, int nclust) 
	{
          double ss_hat = 0.0;
          double wshape, wscl;
	      for(int i=0; i < nobs; i++) {
	           // find cluster observation i belongs to
	           for(int h=0; h < nclust; h++) {
	               if(labels[i] == h+1) {
	                    ss_hat += (rr[i+1] - locations[h])*(rr[i+1] - locations[h]);
	                    break;
	               }
	           }
	       }
	       wshape = sigdf/2 + nobs/2;
	    //   uu = sqrt(ss_hat/nobs);
	       wscl = 1/(ss_hat/2 + (kappa*sigdf)/2);
	       sig[0] = sqrt(1/rgamma(wshape, wscl));
	}
private:
	//prior parameters
	double psi1;
	double psi2;
};
#endif
