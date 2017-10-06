#include "mex.h"
#include "matrix.h"

void tcdf_call_matlab_clean(double *x, mwSignedIndex N, double *nu, double *cdf)
{
    mxArray *rhs[2], *lhs[1];
    double *xp, *nup, *cdfp;
    int i;
    
    rhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    rhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);   // start copying at the double returned by mxGetPr(array_ptr)
     
    xp = mxGetPr(rhs[0]);
    nup = mxGetPr(rhs[1]);

    for( i=0; i<N; i++ ) {
        xp[i] = x[i];
    }
    nup[0] = nu[0]; 
            
    mexCallMATLAB( 1, lhs, 2, rhs, "tcdf" );
    
    mxDestroyArray(rhs[0]);

    cdfp = mxGetPr(lhs[0]);
    for( i=0; i<N; i++ ) {
        cdf[i] = cdfp[i];
    }

    mxDestroyArray(lhs[0]);
}

// http://stackoverflow.com/questions/26619947/avoid-copying-an-array-when-using-mexcallmatlab
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *x, *nu, *cdf;
    mwSignedIndex N;  
     
    x = mxGetPr(prhs[0]);  
    nu = mxGetPr(prhs[1]);  
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    
    /* create the output matrix */  
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    cdf = mxGetPr(plhs[0]);
    
    tcdf_call_matlab_clean(x, N, nu, cdf);
}