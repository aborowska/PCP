function param_trans = transform_param_garch(param, mode)
    param_trans = param;
     
    switch mode            
%         case 'est_opt' % transformation for optimization --> unbounded
        case 'opt' % transformation for optimization --> unbounded
            param_trans(1,2) = log(param(1,2));
            param_trans(1,3) = log(param(1,3)/(1-param(1,3)));
            param_trans(1,4) = log(param(1,4)/(1-param(1,4)));
        case 'back'% if mode == 'back' (transform back)
            param_trans(1,2) = exp(param(1,2));
            param_trans(1,3) = exp(param(1,3))/(1+exp(param(1,3)));
            param_trans(1,4) = exp(param(1,4))/(1+exp(param(1,4)));
            param_trans(1,4) = min(0.998, param_trans(1,4));            
    end
end