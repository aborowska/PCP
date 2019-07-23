function lpdf = sktpdf(x, nu, lambda, a, b, tau)
% LOG KERNEL OF THE SK-T DENSITY

    if (nargin < 6)
        logc = gammaln((nu+1)/2) - gammaln(nu/2) - 0.5*log(pi*(nu-2));
        c = exp(logc);
        a = 4*lambda*c*((nu-2)/(nu-1));
    %     loga = log(a);
    %     loga = log(4) + log(lambda) + logc + log(nu-2) - log(nu-1);
    %     a = exp(loga);
        logb = 0.5*log(1 + 3*lambda^2 - a^2);    
        b = exp(logb);

        tau = - a./b;
    end
%         indicator1 = ((data(t)-mu(t))./sqrt(h(t))<-a./b);
%         indicator2 = ((data(t)-mu(t))./sqrt(h(t))>=-a./b);
    indicator1 = (x < tau);
    indicator2 = (x >= tau);
    
%     likelihoods1 =  - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
%     likelihoods2 =  - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
% %         likelihoods1 = log(b) + log(c) - ((nu+1)./2).*log(1+1./(nu-2).*((b.*indicator1.*((data(t)-mu(t))./sqrt(h(t)))+a)./(1-lamda)).^2);
% %         likelihoods2 = log(b) + log(c) - ((nu+1)./2).*log(1+1./(nu-2).*((b.*indicator2.*((data(t)-mu(t))./sqrt(h(t)))+a)./(1+lamda)).^2);         
% %         likelihoods = - 0.5*log(h(t)) + indicator1.*likelihoods1 + indicator2.*likelihoods2; 
%     lpdf = sum(likelihoods1) + sum(likelihoods2);    

    lpdf = NaN(size(x));
    lpdf(indicator1) = - ((nu+1)./2).*log(1 + (((b.*x(indicator1) + a)./(1-lambda)).^2)./(nu-2));
    lpdf(indicator2) = - ((nu+1)./2).*log(1 + (((b.*x(indicator2) + a)./(1+lambda)).^2)./(nu-2)); 
end
