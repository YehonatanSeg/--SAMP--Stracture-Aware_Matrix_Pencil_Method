% s is assumed to be the sorted singular value of a Signal Hankel matrix,
% i.e., H - U*S*V'. 
function [r, true_idx] = get_rank_by_singular_value_gap(S)
    
    [S,I] = sort(S,'descend');
    ratio = S(1:end-1) ./ S(2:end);
    [~,r] = max(ratio);
    true_idx = I(1:r);

end