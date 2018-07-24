function fvec_lex = matrix_to_vector_lex( fmat_lex )
% Maps a matrix that its values are stored in lexicographical order
% to a vector.

[rows,cols] = size(fmat_lex);

fvec_lex = zeros(rows*cols,1);

ind = 0;
for i = rows:-1:1
    ind = ind + 1;
    for j = 1:cols
        l = cols*(ind-1) + j;
        fvec_lex(l) = fmat_lex(i,j);
    end
end


end

