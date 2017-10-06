function [nwse, var_OLS] = NeweyWest2(d)
    [N1,N2] = size(d);
    N = max(N1,N2);
    d = reshape(d,N,1);    
    resid = d-mean(d);
   
    L = floor( 4*(N/100)^(2/9) ) + 1;
    iota = ones(N,1);
    OneToN = (1:N)';

    mAbs_i_minus_j = abs( iota * OneToN' - OneToN *  iota' );
    Weights = (ones(N,N) - mAbs_i_minus_j*(1/L)) .* (mAbs_i_minus_j <= L);

    OmegaHat = (resid * resid') .* Weights;
    % var_OLS = (N/(N-1)) * (1/N) * (iota'*mOmegaHat*iota) *  (1/N) ;
    var_OLS = (iota'*OmegaHat*iota)/((N-1)*N) ;
    nwse = sqrt(var_OLS)';
end