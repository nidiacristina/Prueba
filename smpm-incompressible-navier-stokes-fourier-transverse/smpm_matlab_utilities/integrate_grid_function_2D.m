function I = integrate_grid_function_2D( u, x, z, n, mx, mz );
% I = integrate_grid_function( u, x, z, n, mx, mz );
%
% Computes the GLL quadrature integral of the grid function u defined on a
% cartesian multi-element grid.
%
% Takes 6 arguments:
%
%   u  - Matrix, sized (mz*n)x(mx*n), of the grid function to be integrated.
%   x  - Vector, length mx*n, of grid coordinates for the 2nd dimension.
%   z  - Vector, length mz*n, of grid coordinates for the 1st dimension.
%   n  - Number of collocation points.
%   mx - Number of domains in the 2nd dimension.
%   mz - Number of domains in the 1st dimension.
%
% Returns 1 value:
%
%   I - Vector, length mx*n, of the column-wise integrated grid functions.
%

%  7 July 2015
%  Sumedh Joshi
%  Cornell University

% Get some quadrature weights.
[junk, w, P] = lglnodes( double( n ) - 1 );

% Get some element widths.
Lx = max( x ) - min( x );
Lz = max( z ) - min( z );
hx = Lx / double( mx );
hz = Lz / double( mz );

% Loop through grid function, mutiplying by quadrature weights.
for ii = 1:mx
    for jj = 1:n
        for kk = 1:mz
            for ll = 1:n
                ndx       = (ii - 1) * n * n * mz + (jj  - 1) * n * mz + (kk - 1) * n + ll;
                u(ndx, :) = w(ll) .* w(jj) .* u(ndx, :) .* (hx / 2) .* (hz / 2);
            end
        end
    end
end

% Sum the kernel.
I = sum( u );
