function M=vandermonde_x(X_ini, X_compu,n)

M=zeros(n,n);
for i=1:n
    for j=1:n
        M(i,j)=lagrange_basis(X_ini,j,X_compu(i));
    end
end




