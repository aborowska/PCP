function f_out = volatility_skt_gas(theta, y, H)
    T = length(y); 
    M = size(theta,1);
     
    lambda = theta(:,1);
    nu = theta(:,2);
    mu = theta(:,3);
    omega = theta(:,4);
    A = theta(:,5);
    B = theta(:,6); 
    
     
    logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
    c = exp(logc);
    a = 4.*lambda.*c.*((nu-2)./(nu-1));
    logb = 0.5.*log(1 + 3.*lambda.^2 - a.^2);    
    b = exp(logb);
    tau = - a./b;
    
    f = log(var(y)); %omega;
     
    if (nargin == 2)
        H = 1;
    end 
    
    if (H > 0)
        f_out = zeros(M,H);    
    end
    
    for jj = 2:T
        h = exp(f);               
        scale = sqrt(h);        
        z = (y(jj-1,1)-mu(:,1))./scale;
        
        ind_tau = 2*(z >= tau(:,1)) - 1;            

        nom = (nu(:,1)+1).*b(:,1).*z.*(b(:,1).*z+a(:,1));
        den = (nu(:,1)-2).*(1+ind_tau.*lambda(:,1)).^2 + (b(:,1).*z+a(:,1)).^2;
        s = 0.5*(nom./den - 1);
        f = omega(:,1) + A(:,1).*s + B(:,1).*f;       % new f at time jj       

        if (jj > T-H)
            f_out(:,jj-(T-H)) = f;
        end
    end
    
    if (H == 0)
        f_out = f;
    end
end