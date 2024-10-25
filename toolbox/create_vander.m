% Creates a vandermonde matrix of size length(v) * dim of the form
% [1 v1 ... (v1)^dim
%  1 v2 ... (v2)^dim
%  . 
%  .
%  1 vm ... (vm)^dim ]

function V = create_vander(v,dim)
    
    % v is expected to be a column vector
    if size(v,2) ~=1
        v = v.';
    end

    V = v.^(0:dim-1);

end
