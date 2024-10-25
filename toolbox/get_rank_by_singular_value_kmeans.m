% s is assumed to be the sorted singular value of a Signal Hankel matrix,
% i.e., H - U*S*V'. 
function [r, true_idx] = get_rank_by_singular_value_kmeans(S)

    idx = kmeans(S,2,'Replicates', 1);

    % the average absolute amplitude of the two subsets

    avg(1) = mean(S( idx == 1 ));
    avg(2) = mean(S( idx == 2 ));
    [~,true_group] = max(avg);

    true_idx = find(idx == true_group);
    r = length(true_idx);
end