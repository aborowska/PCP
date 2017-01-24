function jaco_inv = jacobian_garch(param_trans)
    jaco_inv = diag([1,...
        1/exp(param_trans(1,2)),...
        (1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3))),...
        (1+exp(param_trans(1,4)))*(1+exp(-param_trans(1,4)))]);       
 
end
    