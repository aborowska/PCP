function param_trans = transform_param_garch(param, mode)
    % order: [mu, omega, alpha, beta]
    param_trans = param;
     
    if strcmp(mode,'opt')      % transformation for optimization --> unbounded
        param_trans(1,2) = log(param(1,2));
        param_trans(1,3) = log(param(1,3)/(1-param(1,3)));
        param_trans(1,4) = log(param(1,4)/(1-param(1,4)));
    elseif strcmp(mode,'back') %  transform back 
        param_trans(1,2) = exp(param(1,2));
        param_trans(1,3) = exp(param(1,3))/(1+exp(param(1,3)));
        param_trans(1,4) = exp(param(1,4))/(1+exp(param(1,4)));
        param_trans(1,4) = min(0.998, param_trans(1,4));         
    else 
        error('Unknown transformation mode.')
    end
end