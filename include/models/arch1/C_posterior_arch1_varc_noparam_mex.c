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

void prior_arch1(double *theta, mwSignedIndex N,
        mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if (theta[i+1*N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+3*N] <= 0 ) || (theta[i+3*N] >= 1)) // 0<alpha<1 
        {
            r1[i] = 0;
        }
        if (r1[i] == 1)
        {
            r2[i] = 0;             // unif 
        }
    }
}

/* ********************************************************************* */

void C_posterior_arch1_varc_noparam_mex(double *theta, double *y, double *threshold, double *y_S,
        mwSignedIndex N, mwSignedIndex T, double *d, double *h)
{
    mwSignedIndex *r1;
    double *r2, mu, sigma, gam; 
    double cdf, y_mean, THR;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
     
    prior_arch1(theta, N, r1, r2);

    /* Initialise */ 
    y_mean = 0.0;
    for (j=0; j<T; j++)
    {
        y_mean = y_mean + y[j];                     
    }    
    y_mean = y_mean/T;
    
    if (y_S[0] > 0.0)
    {
        for (i=0; i<N; i++)
        { 
            h[i] = y_S[0]; 
        }    
    }
    else
    {
        for (i=0; i<N; i++)
        {           
            gam = theta[i] - theta[i+2*N];
            h[i] = theta[i+1*N] + theta[i+3*N]*gam*gam/(1-theta[i+3*N]);                          
        }         
    }
    
    /* PDF */
    for (i=0; i<N; i++) 
    {           
        if (r1[i]==1) // compute only for valid draws
        {
            d[i] = r2[i]; // prior
            // the first observation: from the stationary distribution
            THR =  threshold[0]*y_mean;
            if (y[0] < THR)  //pdf
            {
                // stationary distribution for the first observation
                mu = y[0] - theta[i];
                mu = mu*mu;
                sigma = mu/h[i];
                sigma = sigma + log(h[i]); 
                sigma = sigma + log2PI;                                
                d[i] = d[i] - 0.5*sigma; 
            }
            else //cdf
            {
                // stationary distribution for the first observation
                mu = theta[i];
                sigma = sqrt(h[i]);
                normcdf_my_mex(&THR, &mu, &sigma, 1, &cdf);
                d[i] = d[i] + log(1.0-cdf);    
            }
            
            for (j=1; j<T; j++)
            {   
                mu = y[j-1] - theta[i+2*N];
                mu = mu*mu;
                h[i] = theta[i+1*N]*(1.0-theta[i+3*N]) + theta[i+3*N]*mu; 
                
                THR = y[j-1] + y_mean;
                THR = threshold[0]*THR;    
                
                if (y[j] < THR) //pdf
                {
                    mu = y[j] - theta[i];
                    mu = mu*mu;
                    sigma = mu/h[i];
                    sigma = sigma + log(h[i]); 
                    sigma = sigma + log2PI;                                
                    d[i] = d[i] - 0.5*sigma; 
                }
                else //cdf
                {
                    mu = theta[i];                
                    sigma = sqrt(h[i]);
                    normcdf_my_mex(&THR, &mu, &sigma, 1, &cdf);
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
}

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int N, T, *T_out;                              /* size of matrix */
    double *theta, *y, *threshold, *y_S;            /* input*/
    double *d, *h;                                 /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    y_S  = mxGetPr(prhs[3]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
    h = mxGetPr(plhs[2]);
  
    /* call the function */
    C_posterior_arch1_varc_noparam_mex(theta, y, threshold, y_S, N, T, d, h); 
    T_out[0] = T;
}