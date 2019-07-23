function score = C_ScoringRule(dens, cdf, y, threshold)

    [~, T] = size(cdf);
    [D, ~] = size(threshold);
    
    score = zeros(D,T);
    for dd = 1:D
        score(dd,:) = fn_C_ScoringRule(dens, cdf(dd,:), y, threshold(dd,:));
    end
end


function score = fn_C_ScoringRule(dens, cdf, y, threshold)
    % censored, time constant
    ind_C = (y >= threshold);
    ind = (y < threshold);
    score = dens*0; % preallocate
    score(ind) = log(dens(ind));
    score(ind_C) = log(1 - cdf(ind_C));
end