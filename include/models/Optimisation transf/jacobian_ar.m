function jaco_inv = jacobian_ar(param_trans)
    jaco_inv = diag([1,...
        1/exp(param_trans(1,2)),...
        0.5*(1+exp(param_trans(1,3)))*(1+exp(-param_trans(1,3)))]);
end
    