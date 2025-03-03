function [left_mode, lambda, b_MP, b_new] = get_MP_components(Y, rank_type, params)

        Y0 = Y(:,1:end-1);
        Y1 = Y(:,2:end);

          
        r = get_rank(Y0, rank_type, params);
        
        % For lowering noise level in case of full rank
        if r == min(size(Y0,1), size(Y0,2))
            r=r-1;
        end
        
        [U, S, V] = svd(Y0, 'econ');  
        Ur = U(:, 1:r); 
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
    
        Ltilde =  Sr\(Ur')*Y1*Vr;

        [W, D] = eig(Ltilde);
        % left_eigenvectors = (W\(Sr\Ur'));   
        % right_eigenvectors = Vr*W;   
        lambda = diag(D);
    
        left_mode  = Ur*Sr*W;
        right_mode = W \ Vr';       
        
        b_new  = (left_mode(1,:).') .* right_mode(:,1);    

        % the classical MP way
        Vander = create_vander(lambda, sum(size(Y))-1).';
        y = [Y(:,1).',Y(end,2:end)].';
        b_MP = pinv(Vander)*y;

end

function r = get_rank(A, rank_type, params)
    
    [~, S, ~] = svd(A,'econ');
    
    if strcmp(rank_type, 'numeric')
        r = rank(A);
    elseif strcmp(rank_type, 'effective')
        r = round(erank(A));
    elseif strcmp(rank_type, 'cov_effective')
        r = round(erank(cov(A)));
        %     r = round(erank(X1'*X2));
    elseif strcmp(rank_type, 'significant')
        [r, ~] = get_rank_by_significant_digits(diag(S), 2); 
    elseif strcmp(rank_type, 'gap')
        [r, ~] = get_rank_by_singular_value_gap(diag(S));
     elseif strcmp(rank_type, 'kmeans')
         [r, ~] = get_rank_by_singular_value_kmeans(diag(S));

     elseif strcmp(rank_type, 'll_aic')    
         [r, ~,~,~] = get_INFORMATION_CRITERIA_MP_model_order(A, 10, params); 
    elseif strcmp(rank_type, 'll_mdl')    
         [~, r,~,~] = get_INFORMATION_CRITERIA_MP_model_order(A, 10, params);
    elseif strcmp(rank_type, 'll_map')    
         [~, ~,r,~] = get_INFORMATION_CRITERIA_MP_model_order(A, 10, params);
    elseif strcmp(rank_type, 'll_evt')    
         [~, ~,~,r] = get_INFORMATION_CRITERIA_MP_model_order(A, 10, params);
    end


end
