function [draw, accept] = IndMH_mit(mit, kernel,M,BurnIn,GamMat)
% run independent MH algorithm by first drawing from mixture of ts (mit)
% then computing the corresponding importance weights using the density
% kernel
    [draw, lnk] = fn_rmvgt_robust(M+BurnIn, mit, kernel, false);
    lnd = dmvgt(draw, mit, true, GamMat);  
    lnw = lnk - lnd;
    lnw = lnw - max(lnw);
    [ind, a] = fn_MH(lnw);
    draw = draw(ind,:);
    accept = a/(M+BurnIn);
    draw = draw(BurnIn+1:BurnIn+M,:);        
end
