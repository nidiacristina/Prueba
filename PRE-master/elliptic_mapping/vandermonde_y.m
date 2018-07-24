function M=vandermonde_y(Y_ini, Y_compu,n)

M=zeros(n,n);
for i=1:n
    for j=1:n
        M(i,j)=lagrange_basis(Y_ini,i,Y_compu(j));
    end
end
