function r1 = prior_agarch11(theta)
    
%     mu = theta(:,1);
    omega = theta(:,2); % "typo" on purpose: not to confuse with the gamma function
%     mu2 = theta(:,3);
    alpha = theta(:,4);
    beta = theta(:,5);
 
    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);

    r1 = (c1 & c2 & c3);
end
