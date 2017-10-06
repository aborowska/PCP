function score = C_ScoringRule(dens, cdf, y, threshold)

    % censored, time constant
    ind_C = (y >= threshold);
    ind = (y < threshold);
    score = dens*0; % preallocate
    score(ind) = log(dens(ind));
    score(ind_C) = log(1 - cdf(ind_C));
end