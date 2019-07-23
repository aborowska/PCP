function [THETA, mean_theta, median_theta, std_theta, mean_accept] = ...
    RWMH_skt_gas(kernel, prior, mu_init_C, delta, MRW, BurnInRW, plot_on)
    % lambda nu mu omega A B
% prior =  @(xxx) prior_skt_gas(xxx, y, hyper, y_S)

    % MRW = 50000;
    % BurnInRW = 10000;
    %delta = 2.38*std_draw; %    1.7549    0.1104   56.4399    0.0533    0.0551
% clear THETA ACCEPT oldlikhood newlikhood mit_C draw_C accept_C mu_Cb Sigma_Cb
    if isempty(delta)
%         delta = [0.08, 5.5, 0.5, 0.015, 0.045, 0.008]; % var10%
%         delta = [0.08, 5, 0.6, 0.015, 0.04, 0.006]; % var15%
%         delta = [0.1, 6.3, 0.2, 0.02, 0.1, 0.01]; % var20%
%         delta = [0.1, 6.3, 0.2, 0.02, 0.1, 0.01]; % var30%
%         delta = [0.1, 5.0, 0.18, 0.03, 0.13, 0.02]; % var40%
%         delta = [0.14, 15.5, 0.7, 0.03, 0.07, 0.015]; % const5%    
%         delta = [0.08, 6.5, 0.4, 0.02, 0.04, 0.005]; % const10%-20%    
%         delta = [0.08, 6.5, 0.3, 0.02, 0.06, 0.01]; % const30%   
%         delta = [0.085, 4.5, 0.1, 0.02, 0.10, 0.027]; % const40%   
        delta = [0.06, 3, 0.1, 0.02, 0.10, 0.02]; % posterior   
    end


    % mean(ACCEPT)

    theta = mu_init_C; % mu_init_C = mu_C  % mu_C = results_a.mit_C.p*results_a.mit_C.mu 
    oldlikhood = kernel(theta);

    THETA = zeros(MRW,6);
    ACCEPT = zeros(MRW,6);

    tic
    for mm = -BurnInRW:1:MRW
        if mod(mm,1000) == 0
            toc_time = toc;
            fprintf('RWMH iter = %i, time = %6.4f s\n',mm, toc_time);       
        end
        accept = zeros(1,6);

        for ii = 1:6
            % Keep a record of the current theta value being updated
            oldtheta = theta(ii);

            % Propose a new value using a RW with uniform proposal density
            theta(ii) = theta(ii) + delta(ii)*randn; 
            if (isfinite(prior(theta)) && (theta(2)< 50))
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
            THETA(mm,:) = theta;
            ACCEPT(mm,:) = accept;
        end
    end
    
    mean_theta = mean(THETA);
    median_theta = median(THETA);
    std_theta = std(THETA);
    mean_accept = mean(ACCEPT);

    if plot_on
        figure(234)
        for ii = 1:6
            subplot(2,4,ii);
            plot(THETA(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end



        figure(345)
        for ii = 1:6
            subplot(2,4,ii);
            histogram(THETA(:,ii));
%             histogram(draw_C(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end
    end
end