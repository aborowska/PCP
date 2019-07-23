function [mean_theta, median_theta, std_theta, mean_accept, THETA] = ...
    RWMH_tgarch11(kernel, prior, mu_init, delta, MRW, BurnInRW, plot_on)
    %nu mu omega alpha beta
    % kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);

    %delta = 2.38*std_draw; %    1.7549    0.1104   56.4399    0.0533    0.0551
    %    0.7911    0.6469    0.8191    0.0067    0.0096
    % delta = [3.5, 0.3, 100, 0.01, 0.01];
    %     0.6907    0.4462    0.3776    0.1074    0.1143
    % delta = [5.5, 0.33, 100, 0.005, 0.005];
    %     0.5528    0.3892    0.2003    0.2533    0.2724
    % delta = [6.0, 0.33, 90, 0.005, 0.005];
    %     0.5698    0.4240    0.6022    0.0613    0.0625
    % delta = [7.5, 0.33, 100, 0.005, 0.005];
    %    0.5021    0.3826    0.3604    0.1941    0.1953
    % delta = [8.5, 0.33, 100, 0.002, 0.002];
    %    0.5322    0.4494    0.4155    0.1974    0.2035
    % delta = [8.5, 0.33, 100, 0.001, 0.001];
    %     0.4796    0.4102    0.7082    0.1227    0.1223
    if isempty(delta)
        delta = [8.5, 0.1, 20, 0.001, 0.001];
    end
    %     0.4769    0.4008    0.3654    0.3778    0.3866


    % mean(ACCEPT)

    theta = mu_init;
    oldlikhood = kernel(theta);

    % MRW = 50000;
    % BurnInRW = 20000;
    THETA = zeros(MRW,5);
    ACCEPT = zeros(MRW,5);

    tic
    for mm = -BurnInRW:1:MRW
        if mod(mm,1000) == 0
            toc_time = toc;
            fprintf('RWMH iter = %i, time = %6.4f s\n',mm, toc_time);       
        end
        accept = zeros(1,5);

        for ii = 1:5
            % Keep a record of the current theta value being updated
            oldtheta = theta(ii);

            % Propose a new value using a RW with uniform proposal density
            theta(ii) = theta(ii) + delta(ii)*randn; 
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
        for ii = 1:5
            subplot(2,3,ii);
            plot(THETA(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end



        figure(345)
        for ii = 1:5
            subplot(2,3,ii);
            histogram(THETA(:,ii));
%             histogram(draw_C(:,ii));
%             title(params{ii},'Interpreter','Latex');
        end
    end
end