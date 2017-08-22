function h_T = volatility_t_garch(theta, y, y_S)
 
    [N ,~] = size(theta);
    mu = theta(:,2);
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);


    T = size(y,1);
%     ind = 2:T;
    y = y';
    
    h = zeros(N,T); 
    h(:,1) = y_S*ones(N,1);
    
%     h(:,ind) = repmat(omega(:,1),1,T-1) + repmat(alpha(:,1),1,T-1).*(repmat(y(1,ind-1),N,1)-repmat(mu(:,1),1,T-1)).^2;
    for jj = 2:T
        h(:,jj) = omega(:,1).*(1-alpha(:,1)- beta(:,1)) + alpha(:,1).*(y(1,jj-1) - mu(:,1)).^2 + beta(:,1).*h(:,jj-1) ;
    end
    h_T = h(:,T);
end