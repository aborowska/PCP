function r2 = prior_skt_agarch11(theta, hyper)
%     T = length(y);
    M = size(theta,1);
    
    lambda = theta(:,1);
    nu = theta(:,2);
%     mu = theta(:,3);
    omega = theta(:,4);
%     mu2 = theta(:,5);
    alpha = theta(:,6);
    beta = theta(:,7);
 
 
    % uniform prior on alpha and beta on (0,1)
    % with restriction alpha + beta < 1
    % prior is an Nx2 matrix: 
    % 1 col - constraint satisfied?
    % 2 col - prior val at the corresponding point

    c1 = ((alpha >= 1e-6) & (alpha < 1) & (beta >= 1e-6) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 1e-6);
    c4 = (nu > 2);
    c5 = ((lambda > -1) & (lambda < 1));
    
    r1 = (c1 & c2 & c3 & c4 & c5);

    r2 = -Inf*ones(M,1);
    r2(r1==true) = log(0.5) + log(hyper) - hyper*(nu(r1==true) - 2);
    
 end
