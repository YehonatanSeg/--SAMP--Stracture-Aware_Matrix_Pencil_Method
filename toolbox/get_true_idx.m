function true_ind = get_true_idx(epsilon, method)

if strcmp(method, 'kmeans')
     % find the true indices using kmeans
    idx = kmeans(epsilon,2,'Replicates', 1);

    % the average absolute amplitude of the two subsets

    avg(1) = mean(epsilon( idx == 1 ));
    avg(2) = mean(epsilon( idx == 2 ));
    [~,true_group] = max(avg);
    true_ind = find(idx == true_group);

elseif strcmp(method, 'SD')
    [S,I] = sort(epsilon,'descend');
    p=2;
    threshold = 10^(-p);
    below_thresh = find( (S/max(S)) < threshold );
    
    if below_thresh
        r = below_thresh(1)-1;
    else
        % if all the value of S/max(S) >= thresh then all the values are significant
        r = length(S); 
    end
    r = max(r,1);
    true_ind = I(1:r);

elseif strcmp(method, 'GAP')
    [S,I] = sort(epsilon,'descend');
    ratio = S(1:end-1) ./ S(2:end);
    [~,r] = max(ratio);
    true_ind = I(1:r);

end

end