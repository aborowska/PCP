#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;
#define  log2PI 1.837877066409345;

/* ********************************************************************* */

void prior_ar1(double *theta, mwSignedIndex N,
        mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    double beta =  4.632823911110911; // constant for the prior: -log(beta(20,1.5)
    double rho;
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if (theta[i+N] <= 0) // sigma>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+2*N] <= -1 ) || (theta[i+2*N] >= 1)) // -1<rho<1 
        {
            r1[i] = 0;
        }
      
        if (r1[i] == 1)
        {
            r2[i] = -log(theta[i+N]);           // sigma ~ 1/sigma
            rho = (theta[i+2*N]+1.0)/2.0;       // (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);                   
            r2[i] = r2[i] + beta + (20-1)*log(rho) + (1.5-1)*log(1-rho);  
        }
    }
}

/* ********************************************************************* */

void posterior_ar1_mex(double *y, mwSignedIndex N, mwSignedIndex T,
        double *theta, double *d)
{
    mwSignedIndex *r1;
    double *r2, *a1, *P1, mu, sigma; 
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    /* Stationary distibution for the first observation */
    a1 = mxMalloc((N)*sizeof(double));              
    P1 = mxMalloc((N)*sizeof(double));              
     
    prior_ar1(theta, N, r1, r2);

    /* Initialise */     
    for (i=0; i<N; i++)
    {
        a1[i] = theta[i]/(1.0-theta[i+2*N]);      
        P1[i] = (theta[i+N]*theta[i+N])/(1.0 - theta[i+2*N]*theta[i+2*N]);
    }    
         
    /* PDF */
    for (i=0; i<N; i++) 
    {           
        if (r1[i]==1) // compute only for valid draws
        {
            d[i] = r2[i]; // prior
            // the first observation: from the stationary distribution 
            mu = y[0] - a1[i];
            mu = mu*mu;
            sigma = mu/P1[i];
            sigma = sigma + log(P1[i]); 
            sigma = sigma + log2PI;                                
            d[i] = d[i] - 0.5*sigma;
            // Gaussian constants            
            sigma = theta[i+N]*theta[i+N];
            sigma = log(sigma);
            sigma = sigma + log2PI;
            d[i] = d[i] - 0.5*(T-1)*sigma; 
            
            for (j=1; j<T; j++)
            {   
                mu = y[j] - theta[i] - theta[i+2*N]*y[j-1];
                mu = mu*mu;
                mu = mu/(theta[i+N]*theta[i+N]);
                d[i] = d[i] - 0.5*mu;
             }         
        }
        else
        {
                d[i] = -mxGetInf();
        }    
     }
    
    /* Free allocated memory */
    mxFree(r1); 
    mxFree(r2); 
    mxFree(a1); 
    mxFree(P1); 
}

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int N, T, *T_out;                              /* size of matrix */
    double *theta, *y;                             /* input*/
    double *d;                                     /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
  
    /* call the function */
    posterior_ar1_mex(y, N, T, theta, d); 
    T_out[0] = T;
}