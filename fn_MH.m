function [ind, a] = fn_MH(lnw, print_on)
% Run the independence Metropolis-Hastings algorithm on the vector of log weights lnw
% to return the index indicating the corresponding draws  
% a - counter for accepted draws
    if (nargin == 1)
        print_on = false;
    end
    [N,~] = size(lnw);
    a = 0;
    % initialize the chain
    ind = zeros(N,1);
    ind(1,:) = 1;

    ind_old = 1;
    if print_on
        fprintf('\nMH running...')
        % h = waitbar(0,'MH in progress...');
    end
    
    % iterate over the chain
    for ii = 2:N
        u = rand; % draw from uniform
%         e = min(1,w(ii)/w(ind_old)); % min between the ratio and 1
        e = min(1, exp(lnw(ii)-lnw(ind_old))); % min between the ratio and 1
        if (u <= e)
            ind(ii,:) = ii;
            a = a + 1;  % increase the counter for acceptance rate
            ind_old = ii; % move the old position
        else 
            ind(ii,:) = ind_old;
        end 
%         if print_on
%          waitbar(ii/N)
%         end
    end
    if print_on
        fprintf(' done! \n')
        %     close(h) 
    end
end