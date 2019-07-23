function results = PCP_CP_agarch11_estimate_evaluate(kernel, kernel_init, y, draw, ...
    p_bar0, p_bar1, p_bar, T, H, ...
    sampling_opt, cont, options, GamMat)


    y_S = var(y(1:T));    
 
    M = sampling_opt. M;
    BurnIn = sampling_opt.BurnIn;
    mu_init = sampling_opt.mu_init;
    BurnIn_PCP = sampling_opt.BurnIn_PCP;
    thinning = sampling_opt.thinning;
    df = sampling_opt.df;   
    partition = sampling_opt.partition;
    II = sampling_opt.II;



    fprintf('*** Censored Posterior ***\n');

%     cont.mit.CV_tol = 0.3; 
    cont.mit.CV_max = 1.9;
    CV_C = cont.mit.CV_old;
    while (CV_C(end) >= 2)
        try
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            Sigma_C = inv(T*Sigma_C);
            draw_C = rmvt(mu_C,Sigma_C,df,M+BurnIn);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);
        catch
            [mu_C,~,~,~,~,Sigma_C] = fminunc(kernel_init,mu_init,options);
            mit_C = struct('mu',mu_C,'Sigma',reshape(Sigma_C,1,length(mu_C)^2),'df', df, 'p', 1);
            [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            if CV_C(end)>2
                [mit_C, CV_C] = MitISEM_new2(mit_C, kernel, mu_init, cont, GamMat);   
            end
            [draw_C, accept_C] = IndMH_mit(mit_C, kernel,M,BurnIn,GamMat);  
        end 
    end
    h_post_C = volatility_agarch11(draw_C,y,y_S,H);
    y_post_C = randn(M,H).*sqrt(h_post_C);
    y_post_C = bsxfun(@plus,y_post_C,draw_C(:,1));
    y_post_C = sort(y_post_C);

    VaR_05_post_C = y_post_C(p_bar0*M,:); 
    VaR_1_post_C = y_post_C(p_bar1*M,:); 
    VaR_5_post_C = y_post_C(p_bar*M,:); 

    ES_05_post_C = mean(y_post_C(1:p_bar0*M,:)); 
    ES_1_post_C = mean(y_post_C(1:p_bar1*M,:)); 
    ES_5_post_C = mean(y_post_C(1:p_bar*M,:)); 
    

    %% PARTIALLY CENSORED: keep alpha and beta uncensored, then censor mu and sigma
    fprintf('*** Partially Censored Posterior ***\n');
    % mit_C: joint candidate for the joint censored posterior
    draw_short = draw((1:II:M)',:); % thinning - to get hight quality rhos
%     M_short = M/II;
    [draw_PC, a_PC] = sim_cond_mit_MH_outloop(mit_C, draw_short,...
        partition, II, BurnIn_PCP, kernel, GamMat, cont.disp, thinning);
    accept_PC = mean(a_PC);   
    
    r1 = prior_agarch11(draw_PC);
    M_fin = sum(r1);
    draw_PC = draw_PC(r1,:);


    h_post_PC = volatility_agarch11(draw_PC,y,y_S,H);
%     y_post_PC = bsxfun(@times,randn(M_fin,H),sqrt(h_post_PC(T+1:T+H,1))');
    y_post_PC = randn(M_fin,H).*sqrt(h_post_PC);
    y_post_PC = bsxfun(@plus,y_post_PC,draw_PC(:,1));
    y_post_PC = sort(y_post_PC);

    VaR_1_post_PC = y_post_PC(round(p_bar1*M_fin),:); 
    VaR_5_post_PC = y_post_PC(round(p_bar*M_fin),:); 
    VaR_05_post_PC = y_post_PC(round(p_bar0*M_fin),:); 

    ES_1_post_PC = mean(y_post_PC(1:round(p_bar1*M_fin),:)); 
    ES_5_post_PC = mean(y_post_PC(1:round(p_bar*M_fin),:)); 
    ES_05_post_PC = mean(y_post_PC(1:round(p_bar0*M_fin),:));
    
    
    results.mit_C = mit_C;
    results.CV_C = CV_C;
    results.draw_C = draw_C;        
    results.accept_C = accept_C;
    results.VaR_C = [VaR_05_post_C;VaR_1_post_C;VaR_5_post_C];
    results.ES_C = [ES_05_post_C;ES_1_post_C;ES_5_post_C];
    
    results.accept_PC = accept_PC;
    results.draw_PC = draw_PC;        
    results.accept_PC = accept_PC;
    results.VaR_PC = [VaR_05_post_PC; VaR_1_post_PC; VaR_5_post_PC];
    results.ES_PC = [ES_05_post_PC; ES_1_post_PC; ES_5_post_PC];
end