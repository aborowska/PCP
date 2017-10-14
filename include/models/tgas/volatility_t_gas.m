function f_out = volatility_t_gas(theta, y, H)
    T = length(y); 
    M = size(theta,1);
       
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    A = theta(:,4);
    B = theta(:,5);
    
%     rho = (nu-2)./nu;    
    nu_con = (nu+1)./(nu-2);
    A = A.*((nu+3)./nu);
    
    f = omega;
     
    if (nargin == 2)
        H = 1;
    end 
    
    if (H > 0)
        f_out = zeros(M,H);    
    end
    
    for jj = 2:T
        C = 1 + ((y(jj-1,1)-mu).^2)./((nu-2).*f);              
        f = omega + A.*(nu_con.*((y(jj-1,1)-mu).^2)./C - f) + B.*f;                        
%         y(jj,1) = mu + sqrt(rho.*f).*trnd(nu);
        if (jj > T-H)
            f_out(:,jj-(T-H)) = f;
        end
    end
    
    if (H == 0)
        f_out = f;
    end
end