function r2 = prior_t_garch11(theta, hyper)
    M = size(theta,1);
    
    nu = theta(:,1);
    mu = theta(:,2);
    omega = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);

    rho = (nu-2)./nu;
 
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);
    c4 = (nu > 2);
    
    r1 = (c1 & c2 & c3 & c4);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*(nu(r1==true) - 2);
end
