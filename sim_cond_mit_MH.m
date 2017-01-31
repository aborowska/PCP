function [draw_partial, accept] = sim_cond_mit_MH(mit, draw, partition, M, BurnIn, kernel, GamMat) %, thinning)
% mit: a structure for a mix of Student's t distributions
% draw: a draw from the joint candidate, ordered in such a way that first
% colums are for the conditional subset, after which the marginal draws are
% stored
% partition: index where the marginal subset starts
% M: number of draws fromthe conditional distribution (for each marginal
% draw)
 
%     if (nargin = 7)
%         thinning = 1;
%     end
    [N,d] = size(draw);
    d1 = d - partition + 1;
    d2 = partition - 1;
    draw_partial = kron(draw,ones(M,1));
    accept = zeros(N,1);
    
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
    d_marg_tot = dmvgt(draw_marg, mit_marg, true, GamMat);
    d_marg = zeros(N,H);
    for h = 1:H
       d_marg(:,h) = log(mit_marg.p(h)) + dmvgt_mex(draw_marg, mit_marg.mu(h), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(1));
    end
    d_marg = bsxfun(@minus, d_marg, d_marg_tot);
    d_marg = exp(d_marg);
%     p_marg = sum(d_marg,2);  % check whether sum up to 1 (for each row)    

% draw components according to their adjusted probabilities

% % DDD = [0.2,0.35,0.45;0.45,0.35,0.2;0.35,0.2,0.45]
% % [scale, ind_comp] = sort(DDD,2)
% % memb_u = rand(3,1)
% % comp = bsxfun(@gt,scale(:,1:end-1),memb_u)
% % memb = 3-sum(comp,2)
% % true_memb = ind_comp(sub2ind([3,3],(1:3)',memb))

    memb = rand(N,1);
    [scale, ind_comp] = sort(d_marg,2);
%     [scale, ind_comp] = sort(d_marg2,2);
    comp = bsxfun(@gt,scale(:,1:end-1),memb);
    memb = H - sum(comp,2);
    true_memb = ind_comp(sub2ind([N,H],(1:N)',memb));

    for h = 1:H
        ind_h = (true_memb == h);
        n_h = sum(ind_h);
        if (n_h>0)
            a_h = zeros(n_h,1);
            draw_h = randn(n_h*(BurnIn+M),d2);
            W = mit_cond.df(h)./chi2rnd(mit_cond.df(h),n_h*(BurnIn+M),1);
            W = sqrt(W);

            Sigma_h12 = mit.Sigma(h,:);
            Sigma_h12 = reshape(Sigma_h12(ind_mix),d1,d2);
            Sigma_h11 = mit.Sigma(h,:);
            Sigma_h11 = inv(reshape(Sigma_h11(ind_marg),d1,d1));

            mu_h_temp = bsxfun(@minus, draw_marg(ind_h,:), mit_marg.mu(h,:));                
            mu_h = bsxfun(@plus, mu_h_temp*Sigma_h11*Sigma_h12,mit_cond.mu(h,:));
  
            lnd_marg =  dmvgt_mex(draw_marg(ind_h,:), mit_marg.mu(h,:), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(1));
   
            Sigma_h22 = mit.Sigma(h,:);
            Sigma_h22 = reshape(Sigma_h22(ind_cond),d2,d2);

            Sigma_h = (Sigma_h22 - Sigma_h12'*Sigma_h11*Sigma_h12)/mit_cond.df(h); % common factor matrix
            
            for jj = 1:n_h % Construct the conditional candidate given theta_1^i
                ind_jj = (1:(BurnIn+M))' + (jj-1)*(BurnIn+M);
    %        Sigma_h is different for each draw
    %        Sigma_h = Sigma_h*(mit.df(h) + mu_h_temp'*Sigma_h11*mu_h_temp); % different scalar factor
                Sigma_hjj = Sigma_h*(mit_marg.df(h) + mu_h_temp(jj,:)*Sigma_h11*mu_h_temp(jj,:)');
                C = chol(Sigma_hjj);
                draw_hjj = draw_h(ind_jj,:)*C;
                draw_hjj = bsxfun(@times,draw_hjj,W(ind_jj,1));  
                % shift by the common mode
                draw_hjj = bsxfun(@plus,draw_hjj,mu_h(jj,:));                  
                % Run the MH algorithm with common marginal draw
                lnk_h = kernel([draw_hjj, repmat(draw_marg(jj,:),(BurnIn+M),1)]);
                % Conditional logdensity 
                lnd_h = dmvgt_mex(draw_hjj, mu_h(jj,:), Sigma_hjj, mit_cond.df(h), 1, GamMat, double(1));
                % Marginal logdensity
                lnd_h = lnd_h + lnd_marg(jj,1);
                lnw_h = lnk_h - lnd_h;
                lnw_h = lnw_h - max(lnw_h);
                [ind_MH, a_h(jj,1)] = fn_MH(lnw_h); % <<<<< MH Algorithm
                draw_h(ind_jj,:) = draw_hjj(ind_MH,:);              
            end
            accept(ind_h,1) = a_h;  
% %             ind_BurnIn = (BurnIn+M)*kron((0:N-1)',ones(M,1)) + kron(ones(M,1),(BurnIn+1:BurnIn+M)');
%             ind_BurnIn = (0:N-1)';
%             ind_BurnIn = ind_BurnIn(ind_h);
            ind_BurnIn = (0:n_h-1)';
            ind_BurnIn = (BurnIn+M)*kron(ind_BurnIn,ones(M,1)) + kron(ones(n_h,1),(BurnIn+1:BurnIn+M)');
%             draw_BurnIn = draw_h(ind_BurnIn,:);
            ind_h = logical(kron(ind_h,ones(M,1)));     
            draw_partial(ind_h,1:d2) = draw_h(ind_BurnIn,:);
        end
    end
    accept = accept/(M+BurnIn);
    
%% Cheat-sheet: drawing from multivariate Student's t:
% Sigma - a scale matrix, m x m 
% C = chol(Sigma);
% X = randn(N,m);
% X = X*C;
% W = df./chi2rnd(df,N,1);
% X = bsxfun(@times,X,sqrt(W));
% X = bsxfun(@plus,X,mu);
end