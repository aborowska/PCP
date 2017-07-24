function draw_partial = sim_cond_inv_trans(draw, grid, kernel)
% draw: a draw from the joint candidate, ordered in such a way that the first
% colum is for the conditional parameter, after which the marginal draws are
% stored

% grid = 2.1:0.1:20;
% N = 180;
% grid = linspace(-1.5,1.5,180)

    M = size(draw,1);
    N = length(grid);

    U = rand(M,1);
    draw_partial= draw;


    for ii = 1:M
        draw_cond = draw(ii, 2:end);
        draw_cond_grid = [grid, repmat(draw_cond,N,1)];
        lnk_grid = kernel(draw_cond_grid);
        cdf_grid = exp(lnk_grid);
        cdf_grid = cumsum(cdf_grid)/sum(cdf_grid);
        k = sum(cdf_grid < u(ii));
        draw_partial(ii,1) = grid(k) + (grid(k+1)-grid(k))*(u-cdf_grid(k))/(cdf_grid(k+1)-cdf_grid(k));    
    end
end