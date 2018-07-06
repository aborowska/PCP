function dd = predict_ar1_hl(theta, y, VaR_prelim)
    T = length(y);
    eps = theta(:,1);
    draw = theta(:,2:end); % mu, sigma, rho
    y_pred = bsxfun(@times,eps,draw(:,2));
    y_pred = bsxfun(@plus,y_pred,draw(:,1));
    y_pred = y_pred + draw(:,3)*y(T);
    ind_not_hl = (y_pred > VaR_prelim);
%     kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold);
%     dd = kernel(theta(:,2:end));
    dd = - 0.5*(log(2*pi) + (theta(:,1).^2)); % normal kernel for the errors
    dd(ind_not_hl) = -inf;
end
