#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932;

/* ********************************************************************* */
 
void sktlogpdf(double x, double nu, double lambda, 
        double a, double b, double tau,  double *logpdf)         
{
    double c, tmp;
    
    c = (nu+1)/2;
    
    if (x < tau)
    {
        tmp = (b*x + a)/(1-lambda);
    }
    else
    {
        tmp = (b*x + a)/(1+lambda);        
    }
    
    tmp = tmp*tmp;
    tmp = tmp/(nu-2);

    logpdf[0] = - c*log(1 + tmp);
}


/* ********************************************************************* */

void prior_skt_agarch_hyper(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
  /*  ORDER: 
theta[i] <= lambda
theta[i+N] <= nu
theta[i+2*N] <= mu
theta[i+3*N] <= omega
theta[i+4*N] <= mu2
theta[i+5*N] <= alpha
theta[i+6*N] <= beta 
 */     
    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if ((theta[i] <= -1) || (theta[i] >= 1)) // -1 < lambda < 1
        {
            r1[i] = 0;
//             mexPrintf("lambda %6.4f\n",theta[i]);
        }     
        if ((theta[i+N] <= 2) || (theta[i+N] > 50) )//nu>2
        {
            r1[i] = 0;
//             mexPrintf("nu %6.4f\n",theta[i+N]);            
        }   
        if ((theta[i+3*N] <= 0) || (theta[i+3*N] > 100))// omega>0
        {
            r1[i] = 0;
//             mexPrintf("omega %6.4f\n",theta[i+3*N]);        
        }          
        if ((theta[i+5*N] <0 ) || (theta[i+5*N] >= 1)) // 0<=alpha<1 
        {
            r1[i] = 0;
//             mexPrintf("alpha %6.4f\n",theta[i+5*N]);        
            
        }
        if ((theta[i+6*N] < 0) || (theta[i+6*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
//             mexPrintf("beta %6.4f\n",theta[i+6*N]);                    
        }
        if (theta[i+5*N] + theta[i+6*N] >= 1) //alpha+beta<1
        {
            r1[i] = 0;
        }
        
        if (r1[i] == 1)
        {
            r2[i] = log(0.5) + log(hyper[0]) - hyper[0]*(theta[i+N] - 2);
//             mexPrintf("hyper %6.4f\n",hyper[0]);        

        }
    }
}

/* ********************************************************************* */

void posterior_skt_agarch11_mex(double *theta, double *y, double *y_S,
        mwSignedIndex N, mwSignedIndex T, double *hyper, 
        double *d, double *h)
{
    
/* ORDER: 
theta[i] <= lambda
theta[i+N] <= nu
theta[i+2*N] <= mu
theta[i+3*N] <= omega
theta[i+4*N] <= mu2
theta[i+5*N] <= alpha
theta[i+6*N] <= beta 
 */       
    
    mwSignedIndex *r1;
    double *r2, gam, z, scale, tau; /* *h */ 
    double logpdf, tmp, mu;
    double a, b, c;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  

    prior_skt_agarch_hyper(theta, hyper, N, r1, r2);

    /* Initialise */ 
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
            gam = theta[i+2*N] - theta[i+4*N];
//             mexPrintf("gam %6.4f\n",gam);
            h[i] = (theta[i+3*N] + theta[i+5*N]*gam*gam)/(1-theta[i+5*N]-theta[i+6*N]);                          
        }         
    }
    
    /* PDF */
    for (i=0; i<N; i++) 
    {           
        if (r1[i]==1) // compute only for valid draws
        {
//             c = lgamma((theta[i+N]+1)/2) - lgamma(theta[i+N]/2) - 0.5*log(PI*(theta[i+N]-2));
            tmp = (theta[i+N]-2);
            tmp = tmp*PI;
//             c = 1.0;
            c = lgamma((theta[i+N]+1)/2) - lgamma(theta[i+N]/2) - 0.5*log(tmp);
            c  = exp(c);
            a = 4*theta[i]*c*((theta[i+N]-2)/(theta[i+N]-1));        // a(pr) = 4.*lambda(pr).*c(pr).*((nu(pr)-2)./(nu(pr)-1));
            tmp = 1 + 3*theta[i]*theta[i] - a*a;
            b = 0.5*log(tmp);       // logb(pr) = 0.5.*log(1 + 3.*lambda(pr).^2 - a(pr).^2);    
            b = exp(b);
        
            tau = -a/b;
            
            d[i] = r2[i]; // prior
            // the first observation: from the stationary distribution               
            scale = sqrt(h[i]);
            z = (y[0] - theta[i+2*N])/scale; // standardised variable
            sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
            //sktpdf(z, nu(ii,1), lambda(ii,1), a(ii,1), b(ii,1), tau(ii,1))      
//             duvt_garch(y[0], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
            d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c));                  
            for (j=1; j<T; j++)
            {   
                mu = y[j-1] - theta[i+4*N]; // mu2
                mu = mu*mu;
                h[i] = theta[i+3*N]*(1.0-theta[i+5*N]-theta[i+6*N]) +
                        theta[i+5*N]*mu + theta[i+6*N]*h[i];               
                
                scale = sqrt(h[i]);
                z = (y[j] - theta[i+2*N])/scale; //mu
                sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
                d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c));       
//                 rhoh = rho[i]*h[i];                   
//                 duvt_garch(y[j], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
//                 d[i] = d[i] + log(pdf);                  
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
    double *theta, *y, *y_S, *hyper;   /* input*/
    double *d, *h;                                 /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
/*  ORDER: 
theta[i] <= lambda
theta[i+N] <= nu
theta[i+2*N] <= mu
theta[i+3*N] <= omega
theta[i+4*N] <= mu2
theta[i+5*N] <= alpha
theta[i+6*N] <= beta 
*/    
    y = mxGetPr(prhs[1]);
    y_S  = mxGetPr(prhs[2]);
    hyper = mxGetPr(prhs[3]); /* hyperparameter on degrees of freedom*/
  
    
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
    posterior_skt_agarch11_mex(theta, y, y_S, N, T, hyper, d, h); 
    T_out[0] = T;
}
