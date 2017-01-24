function theta = rmvgt2_mex(N,mit)
% Random sampling from mixture of t densitites
% N - number of draws
% mit - list with parameters of mixture of t density
    [H, d] = size(mit.mu); % number of components, dimension of t distribution
    
    % sample membership
    memb = randsample(1:H,N,true,mit.p);
    % prepare normal variates
    RND = randn(N,d);

    theta = rmvgt_mex(RND, memb, mit.mu, mit.Sigma, mit.df);
    for h=1:H
        ind_h = (memb == h);
        n_h = sum(ind_h);
        if (n_h>0)
            mu_h = mu(h,:);
            Sigma_h = Sigma(h,:);
            Sigma_h = reshape(Sigma_h,d,d);
            df_h = df(h);
            draw_h = rmvt(mu_h,Sigma_h,df_h,n_h);
            theta(ind_h,:) = draw_h;
        end
    end
end
