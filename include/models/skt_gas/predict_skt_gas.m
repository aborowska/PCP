function [y_hp, eps_hp, f] = predict_skt_gas(theta, y_T, f_T, hp, eps)
    [N ,~] = size(theta);
    
    lambda = theta(:,1);    
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    A = theta(:,5);
    B = theta(:,6); 
    
    logc = NaN*ones(M,1);
    c = NaN*ones(M,1);
    a = NaN*ones(M,1);
    logb = NaN*ones(M,1);
    b = NaN*ones(M,1);
    tau = NaN*ones(M,1);
     
    logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
    c = exp(logc);
    a = 4.*lambda.*c.*((nu-2)./(nu-1));
    logb = 0.5.*log(1 + 3.*lambda.^2 - a.^2);    
    b = exp(logb);
    
%     tau = - a./b;
    
    if (nargin == 4)
         eps_hp = trnd(repmat(nu,1,hp));
    else %(with given eps)
        eps_hp = eps;
    end
    
    y_hp = zeros(N,hp+1);    
    y_hp(:,1) = y_T.*ones(N,1);

    f = zeros(N,hp+1); 
    f(:,1) = f_T;
    
    for jj = 2:(hp+1)
        h = exp(f(:,jj-1));               
        scale = sqrt(h);        
        z = (y_hp(jj,1)-mu(:,1))./scale;
        
        ind_tau = 2*(z >= tau(:,1)) - 1;            

        nom = (nu(:,1)+1).*b(:,1).*z.*(b(:,1).*z+a(:,1));
        den = (nu(:,1)-2).*(1+ind_tau.*lambda(:,1)).^2 + (b(:,1).*z+a(:,1))/^2;
        s = 0.5*(nom./den - 1);
        f(:,jj) = omega(:,1) + A(:,1).*s + B(:,1).*f(:,jj-1);              
%          C = 1 + ((y_hp(:,jj-1)-mu).^2)./((nu-2).*f(:,jj-1));              
%          f(:,jj) = omega + A.*(nu_con.*((y_hp(:,jj-1)-mu).^2)./C - f(:,jj-1)) + B.*f(:,jj-1);                        
% %          f(:,jj) = omega(:,1) + alpha(:,1).*(y_hp(:,jj-1)-mu(:,1)).^2 + beta(:,1).*h(:,jj-1);       
         y_hp(:,jj) = mu(:,1) + sqrt(exp(f(:,jj))).*eps_hp(:,jj-1);
    end
    y_hp = y_hp(:,2:hp+1);
    f = f(:,2:hp+1);
end