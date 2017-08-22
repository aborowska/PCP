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


void tcdf_mex(double *x,  double *mu, double *Sigma, double *nu, mwSignedIndex N,
         double *val)
/* cdf of the univariate t distribution*/
{
    mwSignedIndex i;
    double z;
    mxArray *normInputs[2];
    mxArray *normOutputs[1];
    normInputs[0] = mxCreateDoubleMatrix(1,1,mxREAL); // a double matrix
    // code to fill in normInputs[0]
    normInputs[1] = mxCreateDoubleMatrix(1,1,mxREAL);


//     z = mxMalloc((1)*sizeof(double));    
    
    for (i=0; i<N; i++)
    { 
        z = (x[i]-mu[0])/(Sigma[0]);

        //     /* print out initial matrix */
        //     mexCallMATLAB(0,NULL,1, &x, "disp");
        //     
        //     /* calculate eigenvalues and eigenvectors */
        //     mexCallMATLAB(2, lhs, 1, &x, "eig");
        normInputs[0] = z;
        normInputs[1] = nu[0]; 
        mexCallMATLAB(1,normOutputs,2,normInputs,"tcdf");
        val[i] = normOutputs[0];
    }
    
//     mxFree(z);
mxDestroyArray(normInputs[1]);
mxDestroyArray(normInputs[0]);
// use normOutputs[0]
mxDestroyArray(normOutputs[0]);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N;                  /* size of matrix */
    double *x, *nu, *mu, *Sigma;           /* input*/
    double *val;                      /* output */
    
    /* Getting the inputs */
    x = mxGetPr(prhs[0]);
    mu = mxGetPr(prhs[1]);
    Sigma = mxGetPr(prhs[2]);
    nu = mxGetPr(prhs[3]);

    N = mxGetM(prhs[0]); /* number of draws */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    
    /* get a pointer to the real data in the output matrix */
    val = mxGetPr(plhs[0]);
 
    /* call the function */
    tcdf_mex(x, mu, Sigma, nu, N, val);     
}