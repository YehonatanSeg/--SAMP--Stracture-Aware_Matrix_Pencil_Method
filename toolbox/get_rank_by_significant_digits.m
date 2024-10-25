% Assumes S is ordered from max to min.
function [r, true_idx] = get_rank_by_significant_digits(S, p)
    threshold = 10^-p;
    
    below_thresh = find( (S/max(S)) < threshold );
    
    if below_thresh
        r = below_thresh(1)-1;
    else
        % if all the value of S/max(S) >= thresh then all the values are significant
        r = length(S); 
    end

    r = max(r,1);
    true_idx = 1:r;
end