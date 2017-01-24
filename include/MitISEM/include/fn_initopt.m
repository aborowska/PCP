function [mu, Sigma, val] = fn_initopt(kernel, mu0, cont) % , options)
    d = length(mu0);

    if ~cont.disp
        options = optimset('Display','off');
        id = 'optim:fminunc:SwitchingMethod';
        warning('off',id);
    else
        options = optimset();
%     options = optimset('Display','iter');
%     options = optimset('TolX', 0.0001, 'Display', 'iter', 'Maxiter', 5000, 'MaxFunEvals', 5000, 'LargeScale', 'off', 'HessUpdate', 'bfgs');
    end
    [mu,val,~,~,~,hessian] = fminunc(kernel, mu0 ,options);
%     x = fminsearch(kernel_init,mu_init, options);
%     x = fminsearch(kernel_init,mu_hl), options);
%     [x,~,~,~,~,hessian] = fminunc(kernel,mu0,options)

    try
        [~, T] = kernel(mu0);
        Sigma = inv(T*hessian);
    catch
        Sigma = inv(hessian);
    end
    Sigma = reshape(Sigma,1,d^2);  
end