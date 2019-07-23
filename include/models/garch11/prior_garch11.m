function r1 = prior_garch11(theta)

%     mu = theta(:,1);
    omega = theta(:,2);
    alpha = theta(:,3);
    beta = theta(:,4);
  

    c1 = ((alpha >= 0) & (alpha < 1) & (beta >= 0) & (beta < 1));
    c2 = (alpha + beta < 1);
    c3 = (omega > 0);

    r1 = (c1 & c2 & c3);

%     r2 = -Inf*ones(M,1);
%     r2(r1==true) = log(2);
%     
%     R = [r1, r2];
end
