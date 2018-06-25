function ADG = global_spectral_differentiation_matrix( zlims, n, mz )
% ADG = global_spectral_differentiation_matrix( d, zlims, n, mz )
%
% Assembles the global spectral differentiation matrix in block diagonal form
% and includes the Lax-Wendroff Flux schemen with the flux values set to one.
%
% It is setup to work with the first derivative matrix only and a
% equidistant grid similar to that of smpm_build_cartesian_grid().
%
% Takes 3 arguments:
%
%   zlims - Vector, ~1x2, of vertical domain limits.
%   n     - Number of collocation points per subdomain.
%   mz    - Number of subdomains.
%
% Returns 1 value:
%
%   ADG - Global spectral differentiation matrix of size (n*mz, n*mz).
%

% Gustavo Rivera 2016

xi = lglnodes( double( n ) - 1 );

% Get the spectral differentiation vectors and reshape them into a matrix.
[d, ~, ~] = derv( n - 1, xi , n - 1 );
d         = reshape( d, n, n);

% Allocate the output matrix.
ADG = zeros( n*mz, n*mz );

% Compute subdomain equidistant spacing.
Lz = (zlims(2) - zlims(1)) / double( mz );

% Construct global matrix.
for ks = 1:mz
    k = (ks - 1) * (n);
    for i = 1:n
        for j = 1:n
            ADG(i+k, j+k) = ADG(i+k, j+k) + d(i, j) * 2 / Lz;
        end
    end
end

% Implement flux on LHS and RHS of each block in global matrix.
for ks = 2:mz
    k           = (ks - 1) * n;
    ADG(k+1, k) = -1.0;
    ADG(k, k+1) =  1.0;
end

% Add and subtract flux from each block.
for ks = 2:mz
    k             = (ks - 1) * n;
    ADG(k, k)     = ADG(k, k) - 1 ;
    ADG(k+1, k+1) = ADG(k+1, k+1) + 1 ;
end

end
