#if !defined(_WIN32)
#define dgetrs dgetrs_
#define dgetrf dgetrf_
#endif

#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "blas.h"
#include "lapack.h"

#define  PI  3.1415926535897932;

void dmvt_mex(double *x, double *mu, double *Sigma, double *df,
        double *GamMat,  mwSignedIndex d,  mwSignedIndex G, mwSignedIndex N, 
        double *dens)
/* density of the multivariate Student's t distribution*/
{
    mwSignedIndex i, j, info;
    mwSignedIndex one = 1;
    mwSignedIndex  *IPIV; /* pointer to ints */
    int ind;
    double *xmu;
    double c0, c1, c2, c, e, tmp, etmp, df5;
    double sqrt_det_Sigma;
    char *cht = "T";
    
    /* variable size arrays */  
    xmu = mxMalloc((d)*sizeof(double));    
    IPIV = mxMalloc((d)*sizeof(mwSignedIndex)); 
       
    /* compute the constants */
    c0 = df[0] + d;
    
    if (c0 <= 100)
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
        c1 = GamMat[G-1];
    }
    
    if (df[0] <= 100)
    {
        df5 = df[0]*50000;
        if ((df5 - floor(df5)) < (floor(df5+1) - df5))
        {
            ind = floor(df5);
        }
        else
        {
            ind = floor(df5 + 1);
        }
        ind = ind -  1;
        c2 = GamMat[ind];
    }
    else
    {
        c2 = GamMat[G-1];
    }      
    
    /* Cholesky decomposition of Sigma */
    /* & changes a value into a pointer*/ 
    dgetrf(&d, &d, Sigma, &d, IPIV, &info); 

    /* determinant of Sigma from LU factorisation*/
    sqrt_det_Sigma = 1;
    for (i=0; i<d; i++)
    {
        sqrt_det_Sigma = sqrt_det_Sigma*Sigma[i+d*i];
    }
    sqrt_det_Sigma = pow(fabs(sqrt_det_Sigma),0.5);

    c = df[0]*PI;
    c = pow(c,0.5*d);
    c2 = c2*c;
    c2 = c2*sqrt_det_Sigma;
    c = c1/c2;
    e = -0.5*c0;    

    for (j=0; j<N; j++)
    {
        for (i=0; i<d; i++)
        {
            xmu[i] = x[j+N*i] - mu[i];
        }
        /* inv(Sigma)*xmu'  */
        /* on exit, xmu is the solution X to Sigma*X = xmu i.e. is equal to xmu*inv(Sigma)*/
         dgetrs(cht, &d, &one, Sigma, &d, IPIV, xmu, &d, &info); /* Sigma is the lu factorisation of Sigma */

        tmp = 0; /* tmp = (x-mu)*inv(Sigma)*(x-mu)' */
        for (i=0; i<d; i++)
        {
            tmp = tmp + xmu[i]*(x[j+N*i] - mu[i]); /* xmu = inv(Sigma)*(x-mu)' */
        } 
        tmp = 1 + tmp/df[0];
        etmp = pow(tmp,e);
        dens[j] = exp(log(c) + log(etmp));
    }

    mxFree(IPIV); 
    mxFree(xmu);
}

void dmvgt_mex(double *x, double *mu, double *Sigma, double *df, double *p,
        double *GamMat, double *L,
        mwSignedIndex d,  mwSignedIndex G, 
        mwSignedIndex N, mwSignedIndex H,  
        double *dens)
/* density of the mixture of the multivariate Studenet's t distributions */        

{
    mwSignedIndex i, j;
    double *mu_h, *Sigma_h, *df_h, *dens_h;    

    df_h = mxMalloc((1)*sizeof(double));    
    mu_h = mxMalloc((d)*sizeof(double));    
    Sigma_h = mxMalloc((d*d)*sizeof(double));    
    dens_h = mxMalloc((N)*sizeof(double));
    
    /* initialitsation */
    for (i=0; i<N; i++)
    {
        dens[i] = 0;
    }
    
     for (i=0; i<H; i++)
     {
         df_h[0] = df[i];
         for (j=0; j<(d*d); j++)
         {
             if (j<d)
             {
                 mu_h[j] = mu[i+H*j];
             }
             Sigma_h[j] = Sigma[i+H*j];
         }
         
         dmvt_mex(x, mu_h, Sigma_h, df_h, GamMat, d, G, N, dens_h);  /* density on the i-th component*/

         
         for (j=0; j<N; j++)
         {
             dens[j] = dens[j] + exp(log(p[i]) + log(dens_h[j]));
         }
     }
    if (L[0]==1)
    {
        for (i=0; i<N; i++)
        {
            dens[i] = log(dens[i]);
        }
    }
        
    mxFree(df_h); mxFree(mu_h); mxFree(Sigma_h); mxFree(dens_h);
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *L; /*log or not*/
    mwSignedIndex G, d, N, H;                          /* size of matrix */
    double *x, *mu, *Sigma, *GamMat, *df, *p;          /* input*/
    double *dens;                                      /* output */
    
    /* Getting the inputs */
    x = mxGetPr(prhs[0]);
    mu = mxGetPr(prhs[1]);
    Sigma = mxGetPr(prhs[2]);
    df = mxGetPr(prhs[3]);
    p = mxGetPr(prhs[4]);
    GamMat = mxGetPr(prhs[5]); 
    L = mxGetPr(prhs[6]); 
    
    N = mxGetM(prhs[0]); /* number of draws */
    d = mxGetN(prhs[1]); /* the dimension of t-distribution*/
    H = mxGetM(prhs[1]); /* number of components*/
    G = mxGetM(prhs[5]); /* the length of GamMat*/

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    
    /* get a pointer to the real data in the output matrix */
    dens = mxGetPr(plhs[0]);
           
    /* call the function */
    dmvgt_mex(x, mu, Sigma, df, p, GamMat, L, d, G, N, H, dens);  
    
}
