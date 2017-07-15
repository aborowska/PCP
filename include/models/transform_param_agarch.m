function param_trans = transform_param_agarch(param, mode)
%     % order: [mu, gamma, omega, alpha, beta]
%     theta  = [mu, omega, mu2, alpha, beta]
    
    param_trans = param;
     
    if strcmp(mode,'opt')      % transformation for optimization --> unbounded
%         param_trans(1,3) = log(param(1,3));
        param_trans(1,2) = log(param(1,2));
        param_trans(1,4) = log(param(1,4)/(1-param(1,4)));
        param_trans(1,5) = log(param(1,5)/(1-param(1,5)));
    elseif strcmp(mode,'back') %  transform back 
%         param_trans(1,3) = exp(param(1,3));
        param_trans(1,2) = exp(param(1,2));
        param_trans(1,4) = exp(param(1,4))/(1+exp(param(1,4)));
        param_trans(1,5) = exp(param(1,5))/(1+exp(param(1,5)));
        param_trans(1,5) = min(0.998, param_trans(1,5));         
    else 
        error('Unknown transformation mode.')
    end
end