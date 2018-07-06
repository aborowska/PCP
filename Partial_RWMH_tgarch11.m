function [mean_theta, median_theta, std_theta, mean_accept, THETA_PCP] = ...
    Partial_RWMH_tgarch11(kernel, prior, partition, THETA_post, delta_pcp, BurnInRW, MRW, plot_on)
    % (nu mu omega) | (alpha beta)
   
    if isempty(delta_pcp)
        delta_pcp = [1.0, 0.2, 15];
    end
    
    theta = THETA_post(1,:);
    oldlikhood = kernel(theta);

    M = size(THETA_post,1);
    
    THETA_PCP = zeros(M*MRW,5);
    ACCEPT_PCP = zeros(M*MRW,(partition-1));

    tic
    for dd = 1:M % for every posterior draw ...
        
        if mod(dd,1000) == 0
            toc_time = toc;
            fprintf('Partial RWMH iter = %i, time = %6.4f s\n',dd, toc_time);       
        end
        
%        ind_bloc = 1:MRW; (MRW+1):(2*MRW); (2*MRW+1):(3*MRW); ...; ((M-1)*MRW+1):(M*MRW)
%               d = 1;     2;               3;                 ...;  M
%        ((dd-1)*MRW+1) : (dd*MRW)
        theta = THETA_post(dd,:);
        ind_bloc = ((dd-1)*MRW+1) : (dd*MRW);
        theta_bloc = zeros(MRW,5);
        accept_bloc = zeros(MRW,(partition-1));
       
        for mm = -BurnInRW:1:MRW % ... run a separate MHRW algorithm
         
            accept = zeros(1,(partition-1));

            for ii = 1:(partition-1)
                % Keep a record of the current theta value being updated
                oldtheta = theta(ii);

                % Propose a new value using a RW with uniform proposal density
                theta(ii) = theta(ii) + delta_pcp(ii)*randn; 
                if (isfinite(prior(theta)) && (theta(3)<100))
                    newlikhood = kernel(theta);

                    % Calculate the log(acceptance probability):
                    % Calculate the new likelihood value for the proposed move:
                    % Calculate the numerator (num) and denominator (den) in turn:

                    % The prior terms are already in the kernel value --> not additng them to the acceptance probability 

                    num = newlikhood;
                    den = oldlikhood;

                    % All other prior terms (for other thetas) cancel in the acceptance probability.
                    % Proposal terms cancel since proposal distribution is symmetric.

                    % Acceptance probability of MH step:
                    A = min(1,exp(num-den));
                else
                    A = 0;
                end
                % To do the accept/reject step of the algorithm        
                % Accept the move with probability A:
                if (rand <= A)  % Accept MEAthe proposed move:
                    % Update the log(likelihood) value:
                    accept(ii) = A;
                    oldlikhood = newlikhood;
                else  % Reject proposed move:
                    % theta stays at current value:
                    theta(ii) = oldtheta;
                end
            end  
            if (mm > 0)
                theta_bloc(mm,:) = theta;
                accept_bloc(mm,:) = accept;
            end
            
        end
        ACCEPT_PCP(ind_bloc,:) = accept_bloc;        
        THETA_PCP(ind_bloc,:) = theta_bloc;
    end
    
    
    mean_theta = mean(THETA_PCP);
    median_theta = median(THETA_PCP);
    std_theta = std(THETA_PCP);
    mean_accept = mean(ACCEPT_PCP);

    if plot_on
        figure(234)
        for ii = 1:5
            subplot(2,3,ii);
            plot(THETA_PCP(:,ii));
        end



        figure(345)
        for ii = 1:5
            subplot(2,3,ii);
            histogram(THETA_PCP(:,ii));
        end
    end
end