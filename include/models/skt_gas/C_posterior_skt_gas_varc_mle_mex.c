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
       
void sktcdf(double *x,  mwSignedIndex sum_C, double nu, double lambda, 
        double a, double b, double tau,  double *cdf)
// sktcdf(z, theta[i+N], theta[i],  a, b, tau, &cdf); 
// tcdf_call_matlab(&z, 1, &theta[i+N], &cdf);
{
    double c, *tmp, *tcdfy;
    int i;
     
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

void prior_skt_gas_hyper(double *theta, double *hyper, 
        mwSignedIndex N, mwSignedIndex *r1, double *r2)
{
    mwSignedIndex i;
    
  /*  NEW ORDER: 
 theta[i] <= lambda
 theta[i+N] <= nu
 theta[i+2*N] <= mu
 theta[i+3*N] <= omega
 theta[i+4*N] <= A
 theta[i+5*N] <= B 
 */     
    /* Variable size arrays */
    
    for (i=0; i<N; i++)
    {
        r1[i] = 1;
        if ((theta[i] >= 1) || (theta[i] <= -1)) // -1<=lambda<1
        {
            r1[i] = 0;
        }
        if (theta[i+3*N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+5*N] < 0) || (theta[i+5*N] >= 1)) // 0<=beta<1
        {
            r1[i] = 0;
        }
        if ((theta[i+4*N] < 0) || (theta[i+4*N] >= 1)) // 0<=A<1 ??
        {
            r1[i] = 0;
        }          
        if (theta[i+N] <= 2) //nu>2
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

void C_posterior_skt_gas_varc_mle_mex(double *theta, double *y, 
        double *mu_mle, double *f_mle, double *threshold, 
        mwSignedIndex N, mwSignedIndex T, 
        double *hyper, double *y_S,   
        double *d)
{
        
/*  NEW ORDER: 
 theta[i] <= lambda
 theta[i+N] <= nu
 theta[i+2*N] <= mu
 theta[i+3*N] <= omega
 theta[i+4*N] <= A
 theta[i+5*N] <= B 
 */     
    
    mwSignedIndex *r1, *cond, sum_C, i2;
    double a, b, c, tau;    
    double *r2, mu, sigma, z, *z2; /* *h */ 
    double s, scale, h, *f, tmp, ind_tau, num, den;    
    double logpdf, *cdf, *THR, pdf, sum_cdf;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    /* cond = 1 for "uncensored" observations */
    cond = mxMalloc((T)*sizeof(mwSignedIndex));       
    THR = mxMalloc((T)*sizeof(double));  
    f = mxMalloc((N)*sizeof(double));  
   
    prior_skt_gas_hyper(theta, hyper, N, r1, r2);

    /* Initialise */ 
 
    sum_C = 0;
    for (j=0; j<T; j++)
    { 
        cond[j] = 0;
                            
        THR[j] = (y[j] - mu_mle[2])/exp(f_mle[j]/2);
        if (THR[j] < threshold[0])
        {
            cond[j] = 1;
        }
        else
        {
            sum_C = sum_C + 1;
        }        
        THR[j] = mu_mle[2] + exp(f_mle[j]/2)*threshold[0];
    }

    z2 = mxMalloc((sum_C)*sizeof(double));
    cdf = mxMalloc((sum_C)*sizeof(double));
        
    
    for (i=0; i<N; i++)
    { 
        f[i] = log(y_S[0]); 
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
            c = lgamma((theta[i+N]+1)/2.0) - lgamma(theta[i+N]/2.0) - 0.5*log(tmp);
            c  = exp(c);
            a = 4*theta[i]*c*((theta[i+N]-2.0)/(theta[i+N]-1));        // a(pr) = 4.*lambda(pr).*c(pr).*((nu(pr)-2)./(nu(pr)-1));
            tmp = 1 + 3*theta[i]*theta[i] - a*a;
            b = 0.5*log(tmp);       // logb(pr) = 0.5.*log(1 + 3.*lambda(pr).^2 - a(pr).^2);    
            b = exp(b); 
        
            tau = -a/b;             
//             mexPrintf("a = %6.4f\n", a);               
//             mexPrintf("b = %6.4f\n", b);               
//             mexPrintf("c = %6.4f\n", c);
            d[i] = r2[i]; // prior
            
            h = exp(f[i]);
            scale = sqrt(h);
            z = (y[0] - theta[i+2*N])/scale;
            if (z >= tau)
            {
                ind_tau = 1;
            }
            else
            {
                ind_tau = -1;
            }            
            
            
            // the first observation
            if (cond[0]==1) //pdf
            {
                sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
//                 mexPrintf("pdf[%i]=%16.14f\n",j,logpdf);
//                 mexPrintf("z[%i]=%16.14f\n",0,z);
                
                d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c));             
            }
            else //cdf
            {
                i2 = i2+1;
                z2[i2] = (THR[0] - theta[i+2*N])/scale;               
            }
            
            for (j=1; j<T; j++)
            {   
                num = (theta[i+N]+1)*b*z*(b*z+a);
                den = (theta[i+N]-2)*(1+ind_tau*theta[i])*(1+ind_tau*theta[i]);
                den = den + (b*z+a)*(b*z+a);
                s = 0.5*(num/den - 1);
                f[i] = theta[i+3*N] + theta[i+4*N]*s + theta[i+5*N]*f[i];
                
                h = exp(f[i]);
                scale = sqrt(h);              
                z = (y[j] - theta[i+2*N])/scale; //mu
                
                if (z >= tau)
                {
                    ind_tau = 1;
                }
                else
                {
                    ind_tau = -1;
                }
                
                
                if (cond[j]==1) //pdf
                {
//                     mexPrintf("z[%i]=%16.14f\n",j+1,z);                    
                    sktlogpdf(z, theta[i+N], theta[i], a, b, tau,  &logpdf);  
//                     mexPrintf("pdf[%i]=%16.14f\n",j+1,logpdf);                    
//                     mexPrintf("z[%i]=%16.14f\n",j+1,z);                    
                    d[i] = d[i] + logpdf - log(scale) + (log(b) + log(c));                    
//                     mexPrintf("d[%i]=%16.14f\n",j+1,d[i]);                      
                }
                else //cdf
                {
                    i2 = i2+1;                    
                    z2[i2] = (THR[j] - theta[i+2*N])/scale;                 
                }                                                   
            }
            
            // print d before cdf
            if (sum_C > 0)
            {
                sktcdf(z2, sum_C, theta[i+N], theta[i], a, b, tau, cdf);                
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
    mxFree(THR);
    
 }

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int N, T, *T_out;                              /* size of matrix */
    double *theta, *y, *y_S,  *mu_mle, *f_mle, *threshold, *hyper;   /* input*/
    double *d;                                 /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);    
    y = mxGetPr(prhs[1]);
    mu_mle = mxGetPr(prhs[2]);        
    f_mle = mxGetPr(prhs[3]);        
    threshold = mxGetPr(prhs[4]);
    hyper = mxGetPr(prhs[5]); /* hyperparameter on degrees of freedom*/
    y_S = mxGetPr(prhs[6]); /* hyperparameter on degrees of freedom*/
  
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
//     plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
//     h = mxGetPr(plhs[2]);
  
    /* call the function */
    C_posterior_skt_gas_varc_mle_mex(theta, y, mu_mle, f_mle, 
            threshold, N, T, hyper, y_S, d); 
    T_out[0] = T;
}
