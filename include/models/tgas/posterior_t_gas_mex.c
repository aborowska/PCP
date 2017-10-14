#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;

/******************************************************************* */

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

/*******************************************************************  */


void prior_t_gas_hyper(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= A
 theta[i+4*N] <= B 
 */     

    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        r2[i] = -mxGetInf();
       
        if (theta[i+2*N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+4*N] < 0) || (theta[i+4*N] >= 1)) // 0<=B<1
        {
            r1[i] = 0;
        }
        if (theta[i] <= 2) //nu>2
        {
            r1[i] = 0;
        }     
        if (r1[i] == 1)
        {
            r2[i] = log(hyper[0]) - hyper[0]*(theta[i+4*N] - 2);
        }
    }
}

/*******************************************************************  */


void posterior_t_gas_mex(double *y, mwSignedIndex N, mwSignedIndex T,
        double *theta, double *hyper, double *GamMat, mwSignedIndex G, double *d)
{
    mwSignedIndex *r1;
    double *r2; 
    double *rho, *A, *nu_con;
    double h, *pdf; 
    double rhoh, tmp;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));            
    rho = mxMalloc((N)*sizeof(double));
    A = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));

    pdf = mxMalloc((1)*sizeof(double));
    
    prior_t_gas_hyper(theta, hyper, N, r1, r2);
 /*  OLD ORDER: 
 theta[i] <= mu
 theta[i+N] <= omega
 theta[i+2*N] <= A
 theta[i+3*N] <= B
 theta[i+4*N] <= nu
 */    
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= A
 theta[i+4*N] <= B 
 */ 

    /* Initialise */
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[i]-2)/theta[i];
        A[i] = theta[i+3*N]*(theta[i]+3)/theta[i];
        nu_con[i] = (theta[i]+1)/(theta[i]-2);
    }
    
    /* PDF */
    for (i=0; i<N; i++) 
    {
        if (r1[i]==1)
        {
            d[i] = r2[i];
            h = theta[i+2*N]/(1-theta[i+4*N]); /* omega/(1-B) */

            for (j=1; j<T; j++)
            {   
                tmp = (y[j-1]-theta[i+N])*(y[j-1]-theta[i+N]);
                tmp = tmp/(h*(theta[i]-2));
                tmp = 1 + tmp;
                tmp = nu_con[i]/tmp;
                h = theta[4*N+i]*h + theta[i+2*N] + A[i]*(tmp*(y[j-1]-theta[i+N])*(y[j-1]-theta[i+N]) - h);   
                rhoh = rho[i]*h;
                duvt_garch(y[j], theta[i+N], rhoh, theta[i], GamMat, G, pdf);
                d[i] = d[i] + log(pdf[0]);
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
    mxFree(rho); 
    mxFree(A); 
    mxFree(nu_con); 
    mxFree(pdf);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSignedIndex N, T, G;                              /* size of matrix */
    double *y, *theta, *hyper;                      /* input*/
    double *GamMat;                                     /* input */
    double *d;                                          /* output */

 /*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= A
 theta[i+4*N] <= B 
 */ 
    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    hyper = mxGetPr(prhs[2]); /* hyperparameter on degrees of freedom*/
    GamMat = mxGetPr(prhs[3]);
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[3]);
    
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    
    /* call the function */
    posterior_t_gas_mex(y, N, T, theta, hyper, GamMat, G, d);
  
}
