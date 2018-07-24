function [ fmat_lex ] = vector_to_matrix_lex( fvec_lex,rows,cols )
% Maps a vector that its values are stored in lexicographical order
% to a matrix (in lexicographical order too).

fmat_lex = zeros(rows,cols);

l = 0;
for i = rows:-1:1
    for j = 1:cols
        l = l + 1;
        fmat_lex(i,j) = fvec_lex(l);
    end
end


end

