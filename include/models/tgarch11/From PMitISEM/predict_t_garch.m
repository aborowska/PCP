function [y_hp, eps_hp, h] = predict_t_garch(theta, y_T, h_T, hp)
    [N ,~] = size(theta);
	
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);

      
    rho = (nu-2)./nu;
    
    eps_hp = trnd(repmat(nu,1,hp));
  
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T.*ones(N,1);

    h = zeros(N,hp+1); 
    h(:,1) = h_T;
    
    for jj = 2:(hp+1)
         h(:,jj) = omega(:,1).*(1-alpha(:,1)-beta(:,1)) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);  
         y_hp(:,jj) = mu(:,1) + sqrt(rho(:,1).*h(:,jj)).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
    h = h(:,2:hp+1);
%     y_hp = predict_t_garch_noS_mex(theta, y_T, h_T, eps_hp);
end