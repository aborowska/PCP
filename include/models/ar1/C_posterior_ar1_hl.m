function dd = C_posterior_ar1_hl(theta, y, threshold)
   kernel = @(xx) C_posterior_ar1_mex(xx, y, threshold);
   dd = kernel(theta(:,2:end));
   dd = dd - 0.5*(log(2*pi) + (theta(:,1).^2)); % normal kernel for the errors
end
