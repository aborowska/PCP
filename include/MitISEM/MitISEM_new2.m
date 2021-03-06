function [mit_new, CV] = MitISEM_new2(kernel_init, kernel, mu_init, cont, GamMat)
    d = size(mu_init,2);
 
    N = cont.mit.N;
    Hmax = cont.mit.Hmax;
%     CV_tol = cont.mit.CV_tol;
%     CV_old = cont.mit.CV_old;
    CV_max = cont.mit.CV_max;
    
    norm = cont.mit.norm;
     
    resampl_on = cont.resmpl_on;
    
%% Step 0: Initialization
    % get/define initial mit density shape scale and degrees of freedom
    % aka: naive proposal density
    if isa(kernel_init, 'function_handle')
        [mu, Sigma] = fn_initopt(kernel_init, mu_init, cont);

        try
            r = fn_testSigma(Sigma);
        catch
            r = 1;
        end

        if (r == 1)
            Sigma = ones(d,d) + 2*diag(ones(d,1));
            Sigma = reshape(Sigma,1,d^2);
            if cont.disp
                display('Initial optimzation FAILED.')    
            end
        else
            if cont.disp
                display('Initial optimzation OK.')
            end
        end
        mit_init.mu = mu;
        mit_init.Sigma = Sigma;
        mit_init.df = cont.mit.dfnc;%1;
        mit_init.p = 1;
        
    elseif isa(kernel_init, 'struct')
        mit_init = kernel_init;
    end
    
    % get draws and IS weights from naive  
    [theta, lnk, ~] = fn_rmvgt_robust(N, mit_init, kernel, resampl_on);
%     display(ind_red);
    lnd = dmvgt(theta, mit_init, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);  
    CV = fn_CVstop2(w);
    if cont.disp
        fprintf('CV = ')        
        fprintf('%6.4f ',CV)
        fprintf('\n')
    end
        
%% Step 1: Adaptation 
    % update scale and location using IS with draws from the naive
    % fixing df 
      
    [mu_adapt, Sigma_adapt] = fn_muSigma(theta, w);
    mit_adapt.mu = mu_adapt;
    mit_adapt.Sigma = Sigma_adapt;
    mit_adapt.df = cont.mit.dfnc;
    mit_adapt.p = 1;
    
    % get draws and IS weights from adapted  
    [theta, lnk, ~] = fn_rmvgt_robust(N, mit_adapt, kernel, resampl_on);
%     display(ind_red);
    lnd = dmvgt(theta, mit_adapt, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);
    CV_new = fn_CVstop2(w);
    CV = [CV, CV_new];
    if cont.disp
        fprintf('CV = ')                
        fprintf('%6.4f ',CV)
        fprintf('\n')
    end    
    
%% Step 2: APPLY ISEM
    % optimize mixture using IS weighted EM and get draws from the new mit
    % optimize mode, scale and df
    cont.df.opt = true;
    mit_old = mit_adapt;
    mit_new = fn_optimt(theta, mit_adapt, w, cont, GamMat);

    % get draws and log kernel evaluation
    [theta, lnk, ~] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);
%     display(ind_red);
    lnd = dmvgt(theta, mit_new, true, GamMat);
    w = fn_ISwgts(lnk, lnd, norm);
	
	% stopping criteria
    CV_new = fn_CVstop2(w);
    CV = [CV, CV_new];
    if cont.disp
        fprintf('CV = ')               
        fprintf('%6.4f ',CV)
        fprintf('\n')
    end    
    
    H = length(mit_new.p);  % number of components

%% Step 3: Iterate on the number of mixture components
    % add more mixture components until convergence
    hstop = false;
    iter = 0;
%     while ((H < Hmax) && (hstop == false))
    while ((iter < cont.mit.iter_max) && (hstop == false))
        iter = iter + 1;
        H = H + 1;
        if cont.disp
            fprintf('H = %d\n',H);
        end
        % select the largest weights and corresponding draws
        ind_w = fn_select(w, cont.mit.ISpc);
        theta_nc = theta(ind_w,:);
        w_nc = w(ind_w);

        % NEW COMPONENT
        % compute new component's mode and scale from IS weights

        [mit_nc.mu, mit_nc.Sigma] = fn_muSigma(theta_nc, w_nc);
        mit_nc.df = cont.mit.dfnc;
        mit_nc.p = cont.mit.pnc;
   
        % COBINE OLD AND NC
        % combine the old mixture mit_new and the new component mit_nc
        mit_old = mit_new;
        mit_new = fn_updateMit(mit_new, mit_nc); 

%%% ??? %%%        
        % DRAW FROM COMBINED
        % get draws and log kernel evaluation from the new mixture mit_new 
        [theta, lnk, ~] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);
%         display(ind_red);

        lnd = dmvgt(theta, mit_new, true, GamMat);
        w = fn_ISwgts(lnk, lnd, norm);
%%%%%%%%%%%
        % UPDATE COMBINED
        % update mode, scale and df  of all mixture components
        mit_new = fn_optimt(theta, mit_new, w, cont, GamMat);
        H = size(mit_new.p,2);

        % DRAW FROM UPDATED
        % get new draws from mit and evaluate new IS weights
        [theta, lnk, ~] = fn_rmvgt_robust(N, mit_new, kernel, resampl_on);
%         display(ind_red);   

        lnd = dmvgt(theta, mit_new, true, GamMat);
        w = fn_ISwgts(lnk, lnd, norm);

        % evaluate convergence
%         CV_old = CV(size(CV,2));
        
        [CV_new, hstop_new] = fn_CVstop2(w, CV_max);
        CV = [CV, CV_new];
        if cont.disp
            fprintf('CV = ')               
            fprintf('%6.4f ',CV)
            fprintf('\n')
        end        
%         if (H > 1)
            hstop = hstop_new;
%         end       
    end
  
        
    if (CV(end) > CV(end-1))
        mit_new = mit_old;
    end
end