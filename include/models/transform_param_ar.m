function xx = transform_param_ar(xx,mode)
    % order: [mu, sigma, rho]
    if strcmp(mode,'opt') % transformation for optimization --> unbounded
        xx(:,2) = log(xx(:,2));
        xx(:,3) = log((1+xx(:,3))./(1-xx(:,3)));
    elseif strcmp(mode,'back') % transform back
        xx(:,2) = exp(xx(:,2));
%         xx(:,3) = logsig(xx(:,3));
        xx(:,3) = -1 + 2./(1+exp(-xx(:,3)));
    else 
        error('Unknown transformation mode.')
    end
end