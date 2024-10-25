% f is assumed to be signal of dimentions channels * time
function H = create_hankel_matrix(f,L)

H = zeros(size(f,2)-L , size(f,1)*(L+1) );

for i = 1:size(f,2)-L
    H(i,:) = f(:,i:i+L);
end


end