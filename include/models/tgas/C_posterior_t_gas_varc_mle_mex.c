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
        if (theta[i+2*N] <= 0) // omega>0
        {
            r1[i] = 0;
        }   
        if ((theta[i+4*N] < 0) || (theta[i+4*N] >= 1)) // 0<=beta<1
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

void C_posterior_t_gas_varc_mle_mex(double *theta, double *y, double *mu_mle, double *threshold, 
        mwSignedIndex N, mwSignedIndex T, 
        double *hyper, double *GamMat, mwSignedIndex G,     
        double *d)
{
    
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= A
 theta[i+4*N] <= B 
 */       
    
    mwSignedIndex *r1, *cond, sum_C, i2;
    double *r2, mu, sigma, *z2; /* *h */ 
    double h, *rho, rhoh;
    double rho_mle, A_mle, nu_con_mle;
    double tmp, *A, *nu_con;
    double *cdf, *THR, h_mle, pdf, sum_cdf;
    mwSignedIndex i, j;
        
    /* Variable size arrays */
    r1 = mxMalloc((N)*sizeof(mwSignedIndex));              
    r2 = mxMalloc((N)*sizeof(double));  
    rho = mxMalloc((N)*sizeof(double));   
    A = mxMalloc((N)*sizeof(double));
    nu_con = mxMalloc((N)*sizeof(double));
    /* cond = 1 for "uncensored" observations */
    cond = mxMalloc((T)*sizeof(mwSignedIndex));       
    THR = mxMalloc((T)*sizeof(double));  
    
    prior_t_gas_hyper(theta, hyper, N, r1, r2);

    /* Initialise */ 
    for (i=0; i<N; i++)
    {
        rho[i] = (theta[i]-2)/theta[i];
        A[i] = theta[i+3*N]*(theta[i]+3)/theta[i];
        nu_con[i] = (theta[i]+1)/(theta[i]-2);
    }

    rho_mle = (mu_mle[0]-2)/mu_mle[0];
    A_mle = mu_mle[3]*(mu_mle[0]+3)/mu_mle[0];
    nu_con_mle = (mu_mle[0]+1)/(mu_mle[0]-2);
        
    sum_C = 0;
    for (j=0; j<T; j++)
    { 
        cond[j] = 0;
        if (j==0)
        {
            h_mle = mu_mle[2]/(1-mu_mle[4]); /* omega/(1-B) */
            rhoh = rho[i]*h;
        }
        else
        {
            tmp = (y[j-1]-mu_mle[1])*(y[j-1]-mu_mle[1]);
            tmp = tmp/(h_mle*(theta[i]-2));
            tmp = 1 + tmp;
            tmp = nu_con_mle/tmp;

            h_mle = mu_mle[4]*h + mu_mle[2] + A_mle*(tmp*(y[j-1]-mu_mle[1])*(y[j-1]-mu_mle[1]) - h_mle);   
            rhoh = rho_mle*h_mle;    
        }                            
        THR[j] = (y[j] - mu_mle[1])/sqrt(rho_mle*h_mle);
        if (THR[j] < threshold[0])
        {
            cond[j] = 1;
        }
        else
        {
            sum_C = sum_C + 1;
        }        
        THR[j] = mu_mle[1] + sqrt(rho_mle*h_mle)*threshold[0];
    }

    z2 = mxMalloc((sum_C)*sizeof(double));
    cdf = mxMalloc((sum_C)*sizeof(double));
    
    /* PDF */
    for (i=0; i<N; i++) 
    {           
        if (r1[i]==1) // compute only for valid draws
        {
            i2 = -1;
            d[i] = r2[i]; // prior
            h = theta[i+2*N]/(1-theta[i+4*N]); /* omega/(1-B) */
            rhoh = rho[i]*h;

            // the first observation: from the stationary distribution
            if (cond[0]==1) //pdf
            {
                // stationary distribution for the first observation
//                 rhoh = rho[i]*h[i];
//                 duvt_garch(y[0], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
//                 d[i] = d[i] + log(pdf);     
                
                duvt_garch(y[0], theta[i+N], rhoh, theta[i], GamMat, G, &pdf);
                d[i] = d[i] + log(pdf);                
            }
            else //cdf
            {
                i2 = i2+1;
                // stationary distribution for the first observation
                z2[i2] = (threshold[0] - theta[i+N])/sqrt(rhoh);               
            }
            
            for (j=1; j<T; j++)
            {   
//                 mu = y[j-1] - theta[i+N];
//                 h[i] = theta[i+2*N]*(1.0-theta[i+3*N]-theta[i+4*N]) + theta[i+3*N]*(mu*mu) + theta[i+4*N]*h[i];
                tmp = (y[j-1]-theta[i+N])*(y[j-1]-theta[i+N]);
                tmp = tmp/(h*(theta[i]-2));
                tmp = 1 + tmp;
                tmp = nu_con[i]/tmp;

                h = theta[4*N+i]*h + theta[i+2*N] + A[i]*(tmp*(y[j-1]-theta[i+N])*(y[j-1]-theta[i+N]) - h);   
                rhoh = rho[i]*h;
    
                if (cond[j]==1) //pdf
                {
//                     rhoh = rho[i]*h[i];                   
//                     duvt_garch(y[j], theta[i+N], rhoh, theta[i], GamMat, G, &pdf); /* MU SIGMA NU */
//                     d[i] = d[i] + log(pdf);  
                    duvt_garch(y[j], theta[i+N], rhoh, theta[i], GamMat, G, &pdf);
                    d[i] = d[i] + log(pdf);                      
                }
                else //cdf
                {
                    i2 = i2+1;
                    z2[i2] = (threshold[0] - theta[i+N])/sqrt(rhoh);                 
                }                                                   
            }
            
            // print d before cdf
//             mexPrintf("d[i] = %6.4f\n", d[i]);               
//             mexPrintf("sum_C = %i\n", sum_C);   
            if (sum_C > 0)
            {
                tcdf_call_matlab(z2, sum_C, &theta[i], cdf);
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
    mxFree(rho); 
    mxFree(z2);
    mxFree(cdf);
    mxFree(A); 
    mxFree(nu_con); 
    mxFree(THR);    
 }

/* ********************************************************************* */

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int N, T, G, *T_out;                              /* size of matrix */
    double *theta, *y,  *mu_mle, *threshold, *hyper;   /* input*/
    double *GamMat;       
//     double *d, *h;                                 /* output */
    double *d;                                 /* output */

    /* Getting the inputs */
    theta = mxGetPr(prhs[0]);
/*  NEW ORDER: 
 theta[i] <= nu
 theta[i+N] <= mu
 theta[i+2*N] <= omega
 theta[i+3*N] <= A
 theta[i+4*N] <= B 
 */ 
    
    y = mxGetPr(prhs[1]);
    mu_mle = mxGetPr(prhs[2]);        
    threshold = mxGetPr(prhs[3]);
    GamMat = mxGetPr(prhs[4]);
    hyper = mxGetPr(prhs[5]); /* hyperparameter on degrees of freedom*/
  
    
    N = mxGetM(prhs[0]); /* no of parameter draws */
    T = mxGetM(prhs[1]); /* no. of observations */
    G = mxGetM(prhs[4]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(N,1,mxREAL); 
    plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
//     plhs[2] = mxCreateDoubleMatrix(N,1,mxREAL); 

    /* get a pointer to the real data in the output matrix */
    d = mxGetPr(plhs[0]);
    T_out = mxGetPr(plhs[1]);
//     h = mxGetPr(plhs[2]);
  
    /* call the function */
    C_posterior_t_gas_varc_mle_mex(theta, y, mu_mle, threshold, N, T, hyper, GamMat, G,  d); 
    T_out[0] = T;
}
