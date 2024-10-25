function [left_mode, lambda, b_MP, b_new] = get_MP_components(X, rank_type, params)

        X1 = X(:,1:end-1);
        X2 = X(:,2:end);

        [U, S, V] = svd(X1, 'econ');    
        r = get_rank(X1, rank_type, params);
        
        % for lowering noise level in case of full rank
        if r == size(S,1)
            r=r-1;
        end

        Ur = U(:, 1:r); 
        Sr = S(1:r, 1:r);
        Vr = V(:, 1:r);
    
        Ltilde =  Sr\(Ur')*X2*Vr;

        [W, D] = eig(Ltilde);
        left_eigenvectors = (W\(Sr\Ur'));   
        right_eigenvectors = Vr*W;   
        lambda = diag(D);
    
        left_mode  = Ur*Sr*W;
        right_mode = W \ Vr';       
        
        b_new  = (left_mode(1,:).') .* right_mode(:,1);    

        % the classical MP way
        Vander = create_vander(lambda, sum(size(X))-1).';
        y = [X(:,1).',X(end,2:end)].';
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
     elseif strcmp(rank_type, 'll_mdl')    
         [~, r] = get_AIC_MDL_MP_model_order(A, params);
      elseif strcmp(rank_type, 'll_aic')    
         [r, ~] = get_AIC_MDL_MP_model_order(A, params);
    end


end
