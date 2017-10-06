function var_nw = NeweyWest(theta,L,weighting,p)
    
    if (nargin < 4) 
        p = 0;
    end
    [N1,N2] = size(theta);
    N = max(N1,N2);
    theta = reshape(theta,N,1);
    theta = PreWhiten(theta,p);
    N = size(theta,1);
 
    if ((nargin == 1) || isempty(L))
        L = floor(0.75*N^(1/3)-1);
    end 
    if ((nargin < 3) || isempty(weighting))
        weighting = 'NW';
    end

    if strcmp(weighting,'P') % parzen kernel
        PR_weights = @(aaa) (1 - 6*aaa.^2 + 6*abs(aaa).^3).*(abs(aaa)<=0.5) ...
                    +  (2*(1-abs(aaa)).^3).*(abs(aaa)>0.5).*(abs(aaa)<=1);
        w = PR_weights((1:L)/L);
    elseif strcmp(weighting,'Q') % quadratic spectral
        QS_weights = @(aaa) 25*(sin(6*pi*aaa/5)./(6*pi*aaa/5) - cos(6*pi*aaa/5))./(12*(pi^2).*(aaa.^2));
        w = QS_weights((1:L)/L);
    else % standard NW weights
        w = (L+1-(1:L))/(L+1);
    end
    
    theta = theta - mean(theta);    
    var_nw = sum(theta.^2,1);
    for ii = 1:L
        gamma =  sum(theta(1+ii:N).*theta(1:N-ii),1);
        var_nw = var_nw + w(1,ii)*(2*gamma); 
    end
    var_nw = var_nw/((N-1)*N);
end


function theta_pw = PreWhiten(theta,p)
    if p > 0
        N = length(theta);
        Y = theta(p+1:N);
        X = zeros(N-p,p);
        for jj = 1:p
            X(:,jj) = theta(p+1-jj:N-jj);
        end
        beta = (X'*X)\X'*Y;
        theta_pw = Y - X*beta;       
    else
        theta_pw = theta;
    end
end