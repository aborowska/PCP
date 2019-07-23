#include "mex.h"
#include "math.h"
#include "matrix.h"

#define  PI  3.1415926535897932

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
       
void sktcdf(double *x,  mwSignedIndex sum_C, double nu, double lambda, 
        double a, double b, double tau,  double *cdf)
// sktcdf(z, theta[i+N], theta[i],  a, b, tau, &cdf); 
// tcdf_call_matlab(&z, 1, &theta[i+N], &cdf);
{
    double c, *tmp, *tcdfy;
    int ind, i;
     
    tmp = mxMalloc((sum_C)*sizeof(double));  
    tcdfy = mxMalloc((sum_C)*sizeof(double));  
   
    c = nu/(nu-2);
    c = sqrt(c);
    
    for (i=0; i<sum_C; i++)
    {
        if (x[i] < tau)
        {
            tmp[i] = (b*x[i] + a)/(1-lambda);
        }
        else
        {
            tmp[i] = (b*x[i] + a)/(1+lambda);        
        }

        tmp[i] = c*tmp[i]; 
    }
    tcdf_call_matlab(tmp, sum_C, &nu, tcdfy);    
    
    for (i=0; i<sum_C; i++)
    {
        if (x[i] < tau)
        {
            cdf[i] = (1-lambda)*tcdfy[i];
        }
        else
        {
            cdf[i] = - lambda + (1+lambda)*tcdfy[i];
        }
    }
    
    mxFree(tmp); 
    mxFree(tcdfy); 

}
/* ********************************************************************* */
        
void sktlogpdf(double x, double nu, double lambda, 
        double a, double b, double tau,  double *logpdf)         
{
    double c, tmp;
    int ind;
    
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
        if ((theta[i+3*N] <= 0) || theta[i+3*N]  > 100)// omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+5*N] <0 ) || (theta[i+5*N] >= 1)) // 0<=alpha<1 
        {
            r1[i] = 0;
        }
        if ((theta[i+6*N] < 0) || (theta[i+6*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
        }
        if (theta[i+5*N] + theta[i+6*N] >= 1) //alpha+beta<1
        {
            r1[i] = 0;
        }
        if ((theta[i+N] <= 2) || theta[i+N] > 50)  //nu>2
        {
            r1[i] = 0;
        }   
        if ((theta[i] <= -1) || (theta[i] >= 1)) // -1 < lambda < 1
        {
            r1[i] = 0;
        }           
        if (r1[i] == 1)
        {
            r2[i] = log(0.5) + log(hyper[0]) - hyper[0]*(theta[i+N] - 2);
        }
    }
}

/* ********************************************************************* */

void C_posterior_skt_agarch11_mex(double *theta, double *y, double *threshold, double *y_S,
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
   
    mwSignedIndex *r1, *cond, sum_C, i2;
    double *r2, gam, z, mu, scale, tau, *z2; 
    double a, b, c, tmp;
    double sum_cdf, *cdf, logpdf;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    cond = mxMalloc((T)*sizeof(mwSignedIndex));
     
    prior_skt_agarch_hyper(theta, hyper, N, r1, r2);

    /* Initialise */     
    sum_C = 0;
    for (j=0; j<T; j++)
    {
        cond[j] = 0;
        if (y[j] < threshold[0])
        {
            cond[j] = 1;
        }
        else
        {
            sum_C = sum_C + 1;
        }                
    }    

    z2 = mxMalloc((sum_C)*sizeof(double));
    cdf = mxMalloc((sum_C)*sizeof(double));

    
//     mexPrintf("sum_cond %6.4f\n",sum_cond);
    
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
            i2 = -1;
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
           
            if (cond[0]==1) //pdf
            {
                // stationary distribution for the first observation
                z = (y[0] - theta[i+2*N])/scale; // standardised variable
                sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
//                 mexPrintf("logpdf[%i]=%6.4f\n",1,logpdf);
                
                d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c)); 
            }
            else //cdf
            {
                i2 = i2+1;                
                // stationary distribution for the first observation
                z = (*threshold - theta[i+2*N])/scale;       
                z2[i2] = z;
//                 sktcdf(z, theta[i+N], theta[i],  a, b, tau, &cdf); 
//                 d[i] = d[i] + log(1.0-cdf);     
            }
            
            
            for (j=1; j<T; j++)
            {   
                mu = y[j-1] - theta[i+4*N]; // mu2
                mu = mu*mu;
                h[i] = theta[i+3*N]*(1.0-theta[i+5*N]-theta[i+6*N]) + theta[i+5*N]*mu + theta[i+6*N]*h[i];               
               
                scale = sqrt(h[i]);              

                if (cond[j]==1) //pdf
                {
                    z = (y[j] - theta[i+2*N])/scale; //mu
                    sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
                    d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c)); 
                }
                else //cdf
                {
                    i2 = i2+1;
                    z = (*threshold - theta[i+2*N])/scale;
                    z2[i2] = z;
//                     sktcdf(z, theta[i+N], theta[i], a, b, tau, &cdf);                         
//                     d[i] = d[i] + log(1.0-cdf);                    
                }                 
            }  
            
            
            if (sum_C > 0)
            {
                sktcdf(z2, sum_C, theta[i+N], theta[i], a, b, tau, cdf);                                         
//                 tcdf_call_matlab(z2, sum_C, &theta[i], cdf);
                sum_cdf = 0;
                for (j=0; j<sum_C; j++)
                {
                    sum_cdf = sum_cdf + cdf[j];
                    d[i] = d[i] + log(1.0- cdf[j]);   
                }            
//                 mexPrintf("sum_cdf = %6.4f\n", sum_cdf);  
//                 mexPrintf("d[i]+cdf = %6.4f\n", d[i]);  
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
    mxFree(z2);
    mxFree(cdf);    
}

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
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
   
    int N, T, *T_out;                               /* size of matrix */
    double *theta, *y, *threshold, *y_S, *hyper;    /* input*/
    double *d, *h;                                  /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
    y = mxGetPr(prhs[1]);
    threshold = mxGetPr(prhs[2]);
    y_S  = mxGetPr(prhs[3]);
    hyper = mxGetPr(prhs[4]); /* hyperparameter on degrees of freedom*/
    
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
    C_posterior_skt_agarch11_mex(theta, y, threshold, y_S, N, T, hyper, d, h); 
    T_out[0] = T;
}