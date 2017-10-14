function draw_partial = sim_cond_inv_trans(draw, grid, kernel)
% draw: a draw from the joint candidate, ordered in such a way that the first
% colum is for the parameter to be drawn conditionally on the marginal parameters, 
% which are stored in the remaining columns

% grid = 2.1:0.1:20;
% N = 180;
% grid = linspace(-1.5,1.5,180)
% kernel = @(xx) C_posterior_t_garch11_2_mex(xx, y(1:T,1), threshold, y_S,  GamMat, hyper);
% kernel = @(xx) C_posterior_t_garch11(xx, y(1:T,1), threshold, y_S, hyper)
 
    M = size(draw,1);
    N = length(grid);

    U = rand(M,1);
    draw_partial = draw;

    tic
    for ii = 1:M
        if (mod(ii,1000) == 0)
            fprintf('Inv. transf. simulation iter = %i\n',ii)
            toc
        end
        draw_marg = draw(ii, 2:end);
        draw_marg_grid = [grid', repmat(draw_marg,N,1)];
        lnk_grid = kernel(draw_marg_grid);
        lnk_grid = lnk_grid - max(lnk_grid); % robustify
        cdf_grid = exp(lnk_grid);
        cdf_grid = cumsum(cdf_grid)/sum(cdf_grid);
        k = max(1,sum(cdf_grid < U(ii)));
        draw_partial(ii,1) = grid(k) + (grid(k+1)-grid(k))*( U(ii)-cdf_grid(k))/(cdf_grid(k+1)-cdf_grid(k));    
    end
end