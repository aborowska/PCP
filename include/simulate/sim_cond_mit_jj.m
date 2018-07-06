function [draw_hjj, lnd_h_cond] = sim_cond_mit_jj(M, thinning, mu_hjj_temp, mu_hjj, Sigma_h11, Sigma_h, df1, df2, GamMat)
    % df2  = mit_cond.df(h)
    % df1 = mit_marg.df(h);
    Mjj = M*thinning;
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
     
    %% Comput conditional log of marginal candidate density 
    lnd_h_cond = dmvgt_mex(draw_hjj, mu_hjj, Sigma_hjj, df2, 1, GamMat, double(1));     
end
     