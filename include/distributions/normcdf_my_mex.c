#include "mex.h"
#include "math.h"


void erf(double *x)
/* error function,  Abramowitz and Stegun (1964) formula 7.1.26 */
/* from http://www.johndcook.com/blog/cpp_erf/ */
{
   // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double t, y;
    double sign = 1.0;

    // Save the sign of x
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


void normcdf_my_mex(double *x,  double *mu, double *Sigma, mwSignedIndex N,
         double *val)
/* cdf of the univariate normal distribution*/
{
    mwSignedIndex i;
    double z;

//     z = mxMalloc((1)*sizeof(double));    
    
    for (i=0; i<N; i++)
    {
//         z[0] = (x[i]-mu[0])/(sqrt(2)*Sigma[0]);
//         erf(z);
//         val[i] = 0.5*(1+z[0]);
        z = (x[i]-mu[0])/(sqrt(2)*Sigma[0]);
        erf(&z);
        val[i] = 0.5*(1+z);
    }
    
//     mxFree(z);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N;                  /* size of matrix */
    double *x, *mu, *Sigma;           /* input*/
    double *val;                      /* output */
    
    /* Getting the inputs */
    x = mxGetPr(prhs[0]);
    mu = mxGetPr(prhs[1]);
    Sigma = mxGetPr(prhs[2]);
     
    N = mxGetM(prhs[0]); /* number of draws */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    
    /* get a pointer to the real data in the output matrix */
    val = mxGetPr(plhs[0]);
 
    /* call the function */
    normcdf_my_mex(x, mu, Sigma, N, val);     
}