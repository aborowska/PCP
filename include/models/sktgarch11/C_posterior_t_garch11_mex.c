#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;

/* ********************************************************************* */

void tcdf_call_matlab(double *x, mwSignedIndex N, double *nu, double *cdf)
{
    mxArray *rhs[2], *lhs[1];
    double *xp, *nup, *cdfp;
    int i;
    
    rhs[0] = mxCreateDoubleMatrix(N,1,mxREAL);
    rhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); // start copying at the double returned by mxGetPr(array_ptr)
     
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

/* ********************************************************************* */

void duvt_garch(double x, double mu, double sigma, double df, double *GamMat, mwSignedIndex G, double *pdf)
{
    double c0, c1, c2, c, e, tmp, etmp, df5;
    int ind;
    
    c0 = df+1;
    
    if (c0 <= 100)
    {
        ind = floor(c0*50000) - 1;
        c1 = GamMat[ind];
    }
    else
    {
        c1 = GamMat[G-1];
    }
    
    if (df <= 100)
    {
        df5 = df*50000;
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

    c = df*PI;
    c = pow(c,0.5);
    c2 = c2*c;
    c2 = c2*pow(sigma,0.5);
    c = c1/c2;
    e = -0.5*c0; 
    
    tmp = (x-mu)*(x-mu)/sigma;
    tmp = 1 + tmp/df;
    etmp = pow(tmp,e); 
    pdf[0] = exp(log(c) + log(etmp));
}


/* ********************************************************************* */

void prior_t_garch_hyper(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
  /*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= alpha
 theta[i+4*N] <= beta 
 */     
    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if (theta[i+2*N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+3*N] <0 ) || (theta[i+3*N] >= 1)) // 0<=alpha<1 
        {
            r1[i] = 0;
        }
        if ((theta[i+4*N] < 0) || (theta[i+4*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
        }
        if (theta[i+3*N] + theta[i+4*N] >= 1) //alpha+beta<1
        {
            r1[i] = 0;
        }
        if (theta[i] <= 2) //nu>2
        {
            r1[i] = 0;
        }     
        if (r1[i] == 1)
        {
            r2[i] = log(0.5) + log(hyper[0]) - hyper[0]*(theta[i] - 2);
        }
    }
}

/* ********************************************************************* */

void C_posterior_t_garch11_mex(double *theta, double *y, double *threshold, double *y_S,
        mwSignedIndex N, mwSignedIndex T, 
        double *hyper, double *GamMat, mwSignedIndex G,     
        double *d, double *h)
{
    
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= alpha
 theta[i+4*N] <= beta 
 */       
    
    mwSignedIndex *r1, *cond;
    double *r2, mu, sigma; /* *h */ 
    double *rho, rhoh;
    double  cdf, pdf;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    rho = mxMalloc((N)*sizeof(double));   
    /* cond = 1 for "uncensored" observations */
    cond = mxMalloc((T)*sizeof(mwSignedIndex));
    /* Stationary distibution for the first observation */
//     h = mxMalloc((N)*sizeof(double));              
     
    prior_t_garch_hyper(theta, hyper, N, r1, r2);

    /* Initialise */ 
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[i]-2)/theta[i];
    }

    for (j=0; j<T; j++)
    {
        cond[j] = 0;
        if (y[j] < threshold[0])
        {
            cond[j] = 1;
        }                
    }    
    
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
            h[i] = theta[i+2*N];                          
        }         
    }
    
    /* PDF */
    for (i=0; i<N; i++) 
    {           
        if (r1[i]==1) // compute only for valid draws
        {
            d[i] = r2[i]; // prior
            // the first observation: from the stationary distribution
            if (cond[0]==1) //pdf
            {
                // stationary distribution for the first observation
//                 mu = y[0] - theta[i+N];
//                 mu = mu*mu;
//                 sigma = mu/h[i];
//                 sigma = sigma + log(h[i]); 
//                 sigma = sigma + log2PI;                                
//                 d[i] = d[i] - 0.5*sigma; 
    
//                 h[i] = theta[2*N+i]*h + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);                  
                rhoh = rho[i]*h[i];
                duvt_garch(y[0], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
                d[i] = d[i] + log(pdf);                  
            }
            else //cdf
            {
                // stationary distribution for the first observation
//                 mu = theta[i];
//                 sigma = sqrt(h[i]);                               
//                 normcdf_my_mex(threshold, &mu, &sigma, 1, &cdf);    
                mu = (threshold[0] - theta[i+N])/sqrt(rho[i]*h[i]);
                tcdf_call_matlab(&mu, 1, &theta[i], &cdf);
                
                d[i] = d[i] + log(1.0-cdf);                   
            }
            
            for (j=1; j<T; j++)
            {   
                mu = y[j-1] - theta[i+N];
                h[i] = theta[i+2*N]*(1.0-theta[i+3*N]-theta[i+4*N]) + theta[i+3*N]*(mu*mu) + theta[i+4*N]*h[i];

                if (cond[j]==1) //pdf
                {
//                     mu = y[j] - theta[i+N];
//                     mu = mu*mu;
//                     sigma = mu/h[i];
//                     sigma = sigma + log(h[i]); 
//                     sigma = sigma + log2PI;                                
//                     d[i] = d[i] - 0.5*sigma; 
                    
//                     h[i] = theta[2*N+i]*h + theta[i] + theta[N+i]*(y[j-1]-theta[3*N+i])*(y[j-1]-theta[3*N+i]);   
                    rhoh = rho[i]*h[i];                   
                    duvt_garch(y[j], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
                    d[i] = d[i] + log(pdf);  
                }
                else //cdf
                {
//                     mu = theta[i];                
//                     sigma = sqrt(h[i]);                                       
//                     normcdf_my_mex(threshold, &mu, &sigma, 1, &cdf);                    
                    mu = (threshold[0] - theta[i+N])/sqrt(rho[i]*h[i]);
                    tcdf_call_matlab(&mu, 1, &theta[i], &cdf);
                    d[i] = d[i] + log(1.0-cdf);                    
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
    mxFree(cond); 
    mxFree(rho); 

//     mxFree(h); 
}

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int N, T, G, *T_out;                              /* size of matrix */
    double *theta, *y, *threshold, *y_S, *hyper;   /* input*/
    double *GamMat;       
    double *d, *h;                                 /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= alpha
 theta[i+4*N] <= beta 
 */    
    y = mxGetPr(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    y_S  = mxGetPr(prhs[3]);
    GamMat = mxGetPr(prhs[4]);
    hyper = mxGetPr(prhs[5]); /* hyperparameter on degrees of freedom*/
  
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[4]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
    h = mxGetPr(plhs[2]);
  
    /* call the function */
    C_posterior_t_garch11_mex(theta, y, threshold, y_S, N, T, hyper, GamMat, G,  d, h); 
    T_out[0] = T;
}
