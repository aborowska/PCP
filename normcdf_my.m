function p = normcdf_my(x,mu,sigma)
%NORMCDF Normal cumulative distribution function (cdf).
%   P = NORMCDF(X,MU,SIGMA) returns the cdf of the normal distribution with
%   mean MU and standard deviation SIGMA, evaluated at the values in X.
%   The size of P is the common size of X, MU and SIGMA.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.
%
%   Default values for MU and SIGMA are 0 and 1, respectively.
%
%   [P,PLO,PUP] = NORMCDF(X,MU,SIGMA,PCOV,ALPHA) produces confidence bounds
%   for P when the input parameters MU and SIGMA are estimates.  PCOV is a
%   2-by-2 matrix containing the covariance matrix of the estimated parameters.
%   ALPHA has a default value of 0.05, and specifies 100*(1-ALPHA)% confidence
%   bounds.  PLO and PUP are arrays of the same size as P containing the lower
%   and upper confidence bounds.
%
%   [...] = NORMCDF(...,'upper') computes the upper tail probability of the 
%   normal distribution. This can be used to compute a right-tailed p-value. 
%   To compute a two-tailed p-value, use 2*NORMCDF(-ABS(X),MU,SIGMA).
%
%   See also ERF, ERFC, NORMFIT, NORMINV, NORMLIKE, NORMPDF, NORMRND, NORMSTAT.

%   References:
%      [1] Abramowitz, M. and Stegun, I.A. (1964) Handbook of Mathematical
%          Functions, Dover, New York, 1046pp., sections 7.1, 26.2.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2013 The MathWorks, Inc.


% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
    z = -(x-mu)./(sqrt(2)*sigma);
    p = 0.5*erfc(z);
end