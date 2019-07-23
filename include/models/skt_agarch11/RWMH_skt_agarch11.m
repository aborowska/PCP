function [mean_theta, median_theta, std_theta, mean_accept, THETA] = ...
    RWMH_skt_agarch11(kernel, prior, mu_init, delta, MRW, BurnInRW, plot_on)
    %nu mu omega alpha beta
    % kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);

    %delta = 2.38*std_draw; %    1.7549    0.1104   56.4399    0.0533    0.0551

    if isempty(delta)
        delta = [0.02, 8.5, 0.1, 10, 0.1, 0.001, 0.001];
    end
    %     0.4769    0.4008    0.3654    0.3778    0.3866


    % mean(ACCEPT)

    theta = mu_init;
    oldlikhood = kernel(theta);

    % MRW = 50000;
    % BurnInRW = 20000;
    THETA = zeros(MRW,7);
    ACCEPT = zeros(MRW,7);

    tic
    for mm = -BurnInRW:1:MRW
        if mod(mm,1000) == 0
            toc_time = toc;
            fprintf('RWMH iter = %i, time = %6.4f s\n',mm, toc_time);       
        end
        accept = zeros(1,7);

        for ii = 1:7
            % Keep a record of the current theta value being updated
            oldtheta = theta(ii);

            % Propose a new value using a RW with uniform proposal density
            theta(ii) = theta(ii) + delta(ii)*randn; 
            if (isfinite(prior(theta)) && (theta(4)<100) && (theta(2)< 50))
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
        for ii = 1:7
            subplot(2,4,ii);
            plot(THETA(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end



        figure(345)
        for ii = 1:7
            subplot(2,4,ii);
            histogram(THETA(:,ii));
%             histogram(draw_C(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end
    end
end