function [draw_partial, accept, lnw] = sim_cond_mit_MH_outloop(mit, draw,...
    partition, M, BurnIn, kernel, GamMat, print_on, thinning)
% mit: a structure for a mix of Student's t distributions
% draw: a draw from the joint candidate, ordered in such a way that first
% colums are for the conditional subset, after which the marginal draws are
% stored
% partition: index where the marginal subset starts
% M: number of draws from the conditional distribution (for each marginal
% draw)
% BurnIn: for each marginal draw
% Thinning (optional): scaling factor (retain every thinning'th draw, if 1 then keep all) 
    if (nargin == 8)
        thinning = 1;
    end
    [N,d] = size(draw);
    d1 = d - partition + 1;
    d2 = partition - 1;
    draw_partial = kron(draw,ones(M,1));
    lnw = zeros(N*M,1);
%     draw_partial(:,1:d2) = NaN*draw_partial(:,1:d2); % useful for
%     debugging
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
    clear draw

% marginal mixture
    mit_marg = mit;
    mit_marg.mu = mit_marg.mu(:,partition:end);
    mit_marg.Sigma = mit_marg.Sigma(:,ind_marg);    
     
% initial conditional mixture
    mit_cond = mit;
    mit_cond.mu = mit_cond.mu(:,1:partition-1);  
    mit_cond.df = mit_cond.df + d1;
    
% adjustd component probabilities:
    d_marg_tot = dmvgt(draw_marg, mit_marg, true, GamMat); % total probability (on the whole mixture)
    d_marg = zeros(N,H);
    for h = 1:H
       d_marg(:,h) = log(mit_marg.p(h)) + dmvgt_mex(draw_marg, mit_marg.mu(h,:), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(1));
    end
    d_marg = bsxfun(@minus, d_marg, d_marg_tot); % normalise by the total probability
    d_marg = exp(d_marg);

% draw components according to their adjusted probabilities
    true_memb = multinomial_sampling(d_marg);

% MAIN SAMPLING LOOP OVER COMPONENTS: 
    for h = 1:H
        ind_h = (true_memb == h);
        n_h = sum(ind_h); % how many draws from the current h'th component
        if (n_h>0)
            a_h = zeros(n_h,1);
            draw_h = randn(n_h*M,d2);
            lnw_h = zeros(n_h*M,1);
            % Common component parameters
            Sigma_h12 = mit.Sigma(h,:);
            Sigma_h12 = reshape(Sigma_h12(ind_mix),d1,d2);
            Sigma_h11 = mit.Sigma(h,:);
            Sigma_h11 = inv(reshape(Sigma_h11(ind_marg),d1,d1));

            mu_h_temp = bsxfun(@minus, draw_marg(ind_h,:), mit_marg.mu(h,:));                
            mu_h = bsxfun(@plus, mu_h_temp*Sigma_h11*Sigma_h12,mit_cond.mu(h,:));
            
            % Log density evaluation on the common component
            lnd_marg =  dmvgt_mex(draw_marg(ind_h,:), mit_marg.mu(h,:), mit_marg.Sigma(h,:), mit_marg.df(h), 1, GamMat, double(1));
   
            Sigma_h22 = mit.Sigma(h,:);
            Sigma_h22 = reshape(Sigma_h22(ind_cond),d2,d2);

            Sigma_h = (Sigma_h22 - Sigma_h12'*Sigma_h11*Sigma_h12)/mit_cond.df(h); % common factor matrix
 
            % for each marginal draw from the current h'th component (according to the sampled
            % probabilities) draw M conditional draws           
            for jj = 1:n_h % Construct the conditional candidate given theta_1^i
                print_on_jj = false;
                % For each marginal draw simulate M conditional draws
                if (print_on && (mod(jj,100)==0))
                    fprintf('Component h=%d, conditional sampling iter = %d. ',h,jj)
                    print_on_jj = true;
                end
                ind_jj = (1:M)' + (jj-1)*M;
                draw_marg_jj = draw_marg(jj,:);
                mu_hjj = mu_h(jj,:);                % mu_h2
                mu_hjj_temp = mu_h_temp(jj,:);      % (x-mu_h1)
                % Sigma_h is different for each draw
                % Sigma_h = Sigma_h*(mit.df(h) + mu_h_temp'*Sigma_h11*mu_h_temp); % different scalar factor
                lnd_marg_jj = lnd_marg(jj,1);
                [draw_h(ind_jj,:), a_h(jj,1), lnw_h(ind_jj,:)] = sim_cond_mit_MH_jj(M, BurnIn, thinning, ...
                    draw_marg_jj, mu_hjj_temp, mu_hjj, Sigma_h11, Sigma_h, mit_marg.df(h), mit_cond.df(h),...
                    kernel, lnd_marg_jj, GamMat, print_on_jj);                
            end  
            accept(ind_h,1) = a_h;               
            ind_h = logical(kron(ind_h,ones(M,1)));     
            draw_partial(ind_h,1:d2) = draw_h;   
            lnw(ind_h,:) = lnw_h;
        end
    end
end

function true_memb = multinomial_sampling(d_marg)
    [N, H] = size(d_marg);
    memb = rand(N,1);
    [scale, ind_comp] = sort(d_marg,2);
    comp = bsxfun(@gt,scale(:,1:end-1),memb);
    memb = H - sum(comp,2);
    true_memb = ind_comp(sub2ind([N,H],(1:N)',memb));
end