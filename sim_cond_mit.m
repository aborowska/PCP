function draw = sim_cond_mit(mit, draw, partition, GamMat)
% mit: a structure for a mix of Student's t distributions
% draw: a draw from the joint candidate, ordered in such a way that first
% colums are for the conditional subset, after which the marginal draws are
% stored
% partition: index where the marginal subset starts

    [N,d] = size(draw);
    d1 = d - partition + 1;
    d2 = partition - 1;
    
    ind = 1:d^2; 
    ind = reshape(ind,d,d);
    ind_marg = ind(partition:end,partition:end);
    ind_marg = reshape(ind_marg,1,d1^2);
    ind_cond = ind(1:d2,1:d2);
    ind_cond = reshape(ind_cond,1,d2^2);    
    ind_mix = ind(partition:end,1:d2);
    ind_mix = reshape(ind_mix,1,d1*d2);  
    
    H = size(mit.mu,1);
    draw_marg = draw(:,partition:end);
%     draw_cond = draw(:,1:partition-1);
%     clear draw

% marginal mixture
    mit_marg = mit;
    mit_marg.mu = mit_marg.mu(:,partition:end);
    mit_marg.Sigma = mit_marg.Sigma(:,ind_marg);    
% initial conditional mixture
    mit_cond = mit;
    mit_cond.mu = mit_cond.mu(:,1:partition-1);  
    mit_cond.df = mit_cond.df + d1;
    
% adjustd component probabilities:
%%     d_marg_tot = dmvgt(draw_marg, mit_marg, true, GamMat);
%     d_marg = zeros(N,H);
%     for h = 1:H
%        d_marg(:,h) = log(mit_marg.p(h)) + dmvgt_mex(draw_marg, mit_marg.mu(h), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(1));
%     end
%     d_marg = bsxfun(@minus, d_marg, d_marg_tot);
% %     p_marg = sum(exp(d_marg),2);   
%%     
    d_marg_tot2 = dmvgt(draw_marg, mit_marg, false, GamMat);
    d_marg2 = zeros(N,H);
    for h = 1:H
       d_marg2(:,h) = mit_marg.p(h)*dmvgt_mex(draw_marg, mit_marg.mu(h), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(0));
    end
    d_marg2 = bsxfun(@rdivide, d_marg2, d_marg_tot2);
%     p_marg2 = sum(d_marg2,2); % checkwhether sum up to 1 (for each row)

% draw components according to their adjusted probabilities

% % DDD = [0.2,0.35,0.45;0.45,0.35,0.2;0.35,0.2,0.45]
% % [scale, ind_comp] = sort(DDD,2)
% % memb_u = rand(3,1)
% % comp = bsxfun(@gt,scale(:,1:end-1),memb_u)
% % memb = 3-sum(comp,2)
% % true_memb = ind_comp(sub2ind([3,3],(1:3)',memb))

    memb = rand(N,1);
    [scale, ind_comp] = sort(d_marg2,2);
    comp = bsxfun(@gt,scale(:,1:end-1),memb);
    memb = H - sum(comp,2);
    true_memb = ind_comp(sub2ind([N,H],(1:N)',memb));

    for h = 1:H
        ind_h = (true_memb == h);
        n_h = sum(ind_h);
        if (n_h>0)
            draw_h = randn(n_h,d2);
            W = mit_cond.df(h)./chi2rnd(mit_cond.df(h),n_h,1);
            W = sqrt(W);

            Sigma_h12 = mit.Sigma(h,:);
            Sigma_h12 = reshape(Sigma_h12(ind_mix),d2,d1);
            Sigma_h11 = mit.Sigma(h,:);
            Sigma_h11 = inv(reshape(Sigma_h11(ind_marg),d1,d1));

            mu_h_temp = bsxfun(@minus, draw_marg, mit_marg.mu(h,:));                
            mu_h = bsxfun(@plus, mu_h_temp*Sigma_h11*Sigma_h12',mit_cond.mu(h,:));

            Sigma_h22 = mit.Sigma(h,:);
            Sigma_h22 = reshape(Sigma_h22(ind_cond),d2,d2);

            Sigma_h = (Sigma_h22 - Sigma_h12*Sigma_h11*Sigma_h12')/mit_cond.df(h); % common factor matrix
            for jj = 1:n_h
    %        draw_h = mvtrnd(Sigma_h,df_h,n_h);
    %        draw_h = rmvt(mu_h,Sigma_h,df_h,n_h);

    %        Sigma_h is different for each draw
    %        Sigma_h = Sigma_h*(mit.df(h) + mu_h_temp'*Sigma_h11*mu_h_temp); % different scalar factor
                    Sigma_hjj = Sigma_h*(mit.df(h) + mu_h_temp(jj,:)*Sigma_h11*mu_h_temp(jj,:)');
                    C = chol(Sigma_hjj);
                    draw_h(jj,:) = draw_h(jj,:)*C;
                    draw_h(jj,:) = W(jj,1)*draw_h(jj,:);                
            end
    %             draw_h = draw_h + repmat(mu_h,n_h,1);
            draw_h = draw_h + mu_h;
            draw(ind_h,1:d2) = draw_h;
        end
    end

% drawing from multivariate Student's t:
% C = chol(Sigma);
% X = randn(N,m);
% X = X*C;
% W = df./chi2rnd(df,N,1);
% X = bsxfun(@times,X,sqrt(W));
% X = bsxfun(@plus,X,mu);
end