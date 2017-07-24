#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;
#define  log2PI 1.837877066409345;

/* ********************************************************************* */

void myerf(double *x)
/* error function,  Abramowitz and Stegun (1964) formula 7.1.26 */
/* from http://www.johndcook.com/blog/cpp_erf/ */
{
    /*  constants */
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double t, y;
    double sign = 1.0;

    /*  Save the sign of x */
    if (x[0] < 0)
    {
        sign = -1.0;
    }
    x[0] = fabs(x[0]);

    // A&S formula 7.1.26
    t = 1.0/(1.0 + p*x[0]);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x[0]*x[0]);
    x[0] = sign*y;
}

/* ********************************************************************* */

void normcdf_my_mex(double *x,  double *mu, double *Sigma, mwSignedIndex N,
         double *val)
/* cdf of the univariate normal distribution*/
{
    mwSignedIndex i;
    double z;
   
    for (i=0; i<N; i++)
    {
        z = (x[i]-mu[0])/(sqrt(2)*Sigma[0]);
        myerf(&z);
        val[i] = 0.5*(1+z);
    }
}

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
//             rho = (theta[i+2*N]+1.0)/2.0;       // (phi+1)/2 ~ betapdf((phi+1)/2, 20, 1.5);                   
//             r2[i] = r2[i] + beta + (20-1)*log(rho) + (1.5-1)*log(1-rho);  
        }
    }
}

/* ********************************************************************* */

void C_posterior_ar1_varc_mex(double *y, mwSignedIndex N, mwSignedIndex T,
        double *threshold, double *theta, double *d)
{
    mwSignedIndex *r1;
    double *r2, *a1, *P1, mu, sigma; 
    double THR, cdf;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    /* Stationary distibution for the first observation */
    a1 = mxMalloc((N)*sizeof(double));              
    P1 = mxMalloc((N)*sizeof(double));              
     
    prior_ar1(theta, N, r1, r2);

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
            THR =  threshold[0]*a1[i];             
            if  (y[0] < THR) //pdf // below the (un)condiional mean
            {
                // stationary distribution for the first observation
                mu = y[0] - a1[i];
                mu = mu*mu;
                sigma = mu/P1[i];
                sigma = sigma + log(P1[i]); 
                sigma = sigma + log2PI;                                
                d[i] = d[i] - 0.5*sigma; 
            }
            else //cdf
            {
                // stationary distribution for the first observation
                sigma = sqrt(P1[i]);
//                 normcdf_my_mex(threshold, &a1[i], &sigma, 1, &cdf);
                normcdf_my_mex(&THR, &a1[i], &sigma, 1, &cdf);
                d[i] = d[i] + log(1.0-cdf);    
            }
            
            for (j=1; j<T; j++)
            {   
                mu = theta[i] + theta[i+2*N]*y[j-1];
                THR =  threshold[0]*mu; 
                if (y[j] < THR) //pdf // below the condiional mean
                {
                    mu = y[j] - mu;
                    mu = mu*mu;
                    sigma = theta[i+N]*theta[i+N];
                    mu = mu/sigma;                    
                    // Gaussian constants
                    sigma = log(sigma);
                    sigma = sigma + log2PI;
                    d[i] = d[i] - 0.5*(mu + sigma);   
                }
                else //cdf
                {
// //                     normcdf_my_mex(threshold, &mu, &theta[i+N], 1, &cdf);
                    normcdf_my_mex(&THR, &mu, &theta[i+N], 1, &cdf);
                    d[i] = d[i] + log(1-cdf);                    
                }                
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
    double *theta, *y, *threshold;                 /* input*/
    double *d;                                     /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
  
    /* call the function */
    C_posterior_ar1_varc_mex(y, N, T, threshold, theta, d); 
    T_out[0] = T;
}
