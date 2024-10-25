% Calculate the effective rank according to the paper :
% "THE EFFECTIVE RANK: A MEASURE OF EFFECTIVE DIMENSIONALITY"
function r  = erank(A)
    [~,S,~] = svd(A,"econ");
    
    S_norm = sum(abs(diag(S)));
    p = diag(S)/S_norm;
    
    H = -sum(p.*log(p));
    r = exp(H);


end