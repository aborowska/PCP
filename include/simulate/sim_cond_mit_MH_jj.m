function [draw_hjj, a_jj, lnw_h] = sim_cond_mit_MH_jj(M, BurnIn, thinning, draw_marg_jj, mu_hjj_temp, mu_hjj, Sigma_h11, Sigma_h, ...
    df1, df2, kernel, lnd_marg_jj, GamMat, print_on)
    % df2  = mit_cond.df(h)
    % df1 = mit_marg.df(h);
    Mjj = BurnIn + M*thinning;
    d2 = size(Sigma_h,2); 
    
    %% Build individual Sigma_hjj (mu_hjj_temp = x_jj - mu_hjj)
    % mu_h_temp = bsxfun(@minus, draw_marg(ind_h,:), mit_marg.mu(h,:));                
    Sigma_hjj = Sigma_h*(df1 + mu_hjj_temp*Sigma_h11*mu_hjj_temp');  
    
    %% Generate from the Student's t distribution with mu_hjj, Sigma_hjj, df2
%     draw_h = randn(n_h*(BurnIn+M),d2);
    draw_hjj = randn(Mjj,d2);
%     W = mit_cond.df(h)./chi2rnd(mit_cond.df(h),n_h*(BurnIn+M),1);
    W = df2./chi2rnd(df2,Mjj,1);
    W = sqrt(W);  
    C = chol(Sigma_hjj);
    draw_hjj = draw_hjj*C;
    draw_hjj = bsxfun(@times,draw_hjj,W);  
    % shift by the common mode
    draw_hjj = bsxfun(@plus,draw_hjj,mu_hjj);                   
    
    %% Run the MH algorithm with common marginal draw
    lnk_h = kernel([draw_hjj, repmat(draw_marg_jj,Mjj,1)]);
    % Conditional logdensity 
    lnd_h = dmvgt_mex(draw_hjj, mu_hjj, Sigma_hjj, df2, 1, GamMat, double(1));   
    % Marginal logdensity
    lnd_h = lnd_h + lnd_marg_jj;
    lnw_h = lnk_h - lnd_h;
    lnw_h = lnw_h - max(lnw_h);
    % MH Algorithm
    [ind_MH, a_jj] = fn_MH(lnw_h, print_on);
    draw_hjj = draw_hjj(ind_MH,:);           
    draw_hjj = draw_hjj(BurnIn+1:Mjj,:);
    draw_hjj = draw_hjj(1:thinning:(Mjj-BurnIn),:);
    a_jj = a_jj/Mjj;
    
    
    lnw_h = lnw_h(ind_MH,:);           
    lnw_h = lnw_h(BurnIn+1:Mjj,:);
    lnw_h = lnw_h(1:thinning:(Mjj-BurnIn),:);
end
     