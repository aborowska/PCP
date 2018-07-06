function [Theta, Accept] = fn_MH_RW(theta_init, kernel, prior, M, BurnIn, delta)
% A normal random-walk algorithm 
 
    D = size(theta_init,2); 
    Theta = zeros(M,D);
    Accept = zeros(M,D);
    theta = theta_init;
    oldlnk = kernel(theta);

    % Cycle through each element of theta in turn 
    % and propose to update using random walk MH with normal proposal density
    for t = -BurnIn:M
        accept = zeros(1,D);

        for ii = 1:D
            % Keep a record of the current theta value being updated
            oldtheta = theta(ii);

            % Propose a new value for the parameter theta using a RW with normal proposal density
            % theta(ii) = runif(1, theta(ii)-delta(ii), theta(ii)+delta(ii));
            theta(ii) = theta(ii) + delta(ii)*randn; 
            if isfinite(prior(theta))
                % Calculate the log(acceptance probability):
                % Calculate the new lnk value for the proposed move:
                % Calculate the numerator (num) and denominator (den) in turn:
                newlnk = kernel(theta);

                % % For regression coefficients add in prior (Normal) terms to the acceptance probability
                % num = newlikhood - 0.5*(((theta(ii)-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
                % den = oldlikhood - 0.5*(((oldtheta-prior.T_mu(ii))^2)/prior.T_sigma2(ii));
                %   % All other prior terms (for other thetas) cancel in the acceptance probability.
                % Proposal terms cancel since proposal distribution is symmetric.
                num = newlnk;
                den = oldlnk;
                % Acceptance probability of MH step:
                A = min(1,exp(num-den));
            else
                A = 0;
            end
            % To do the accept/reject step of the algorithm        
            % Accept the move with probability A:
            if (rand <= A)  % Accept the proposed move:
                % Update the log(likelihood) value:
                oldlnk = newlnk;     
                %  accept = accept+1;
                accept(ii) = A;
            else  % Reject proposed move:
                % theta stays at current value:
                theta(ii) = oldtheta;
            end           
        end  

        if (t > 0)
            Theta(t,:) = theta;
            Accept(t,:) = accept;
        end
    end
end