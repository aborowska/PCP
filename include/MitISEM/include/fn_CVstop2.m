function [CV_new, hstop] = fn_CVstop2(w, CV_max)
% compute new CV from IS weights
% and indicator to finalize the number of mixture components in MitISEM
%     if std(w) == 0 
%         error('IS weights w are constnat, try increasing number of draws N.')
%     end

    CV_new = std(w)/mean(w);
  
    if ((nargin > 1) && (nargout > 1))
        hstop = (CV_new <= CV_max);
    end
end