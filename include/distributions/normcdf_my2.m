function p = normcdf_my2(x,mu,sigma)
    z = (x - mu) / (sigma * sqrt(2));
    p = 0.5 * (1 + erf_my(z)); 
end