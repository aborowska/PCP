function p = tcdf_my(x,mu,sigma,v)
%TCDF   Student's T cumulative distribution function (cdf).
%   P = TCDF(X,V) computes the cdf for Student's T distribution
%   with V degrees of freedom, at the values in X.

%   References:
%      [1] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.7.
%      [2] L. Devroye, "Non-Uniform Random Variate Generation",
%      Springer-Verlag, 1986
%      [3] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, Section 10.3, pages 144-146.

%   Copyright 1993-2014 The MathWorks, Inc.

    z = (x-mu)./sigma;
    normcutoff = 1e7;
    [N, d] = size(z);

    % Initialize p
    p = -Inf*ones(N,d);
    ind_not_median = (z~=0);
    p(~ind_not_median) = 0.5; % Make the result exact for the median.

    if (sum(ind_not_median) > 0)
        if(v == 1) % Cauchy
            % Special case for Cauchy distribution.  See Devroye pages 29 and 450.
            % Note that instead of the usual Cauchy formula (atan x)/pi + 0.5, 
            % we use acot(-x), which is equivalent and avoids roundoff error.
            p = xpos(cauchy) + acot(-z)/pi; 
        elseif (v > normcutoff) % normal
            p = normcdf_my(z,0 ,1);
        else
            % General case: first compute F(-|x|) < .5, the lower tail.
            % See Abramowitz and Stegun, formulas and 26.7.1/26.5.27 and 26.5.2
            xsq = x.^2;
            % For small v, form v/(v+x^2) to maintain precision
            t = (v < xsq);
            if any(t(:))
                p(t) = betainc(v(t) ./ (v(t) + xsq(t)), v(t)/2, 0.5, 'lower') / 2;
            end
            % For large v, form x^2/(v+x^2) to maintain precision
            t = (v >= xsq);
            if any(t(:))
                p(t) = betainc(xsq(t) ./ (v(t) + xsq(t)), 0.5, v(t)/2, 'upper') / 2;
            end
            % For x > 0, F(x) = 1 - F(-|x|).
            xpos = (x > 0);
            p(xpos) = 1 - p(xpos); % p < .5, cancellation not a problem
        end
    end
end
