function I = integrate_grid_function_1D( u, x, n, mx );
% I = integrate_grid_function( u, x, n, mx );
%
% Computes the GLL quadrature integral of the grid function u defined on a
% cartesian multi-element grid in 1D.
%
%  Takes 4 arguments:
%
%    u  - Vector, length mx*n, of the grid function to be integrated.
%    x  - Vector, length mx*n, of grid coordinates.
%    n  - Number of collocation points.
%    mx - Number of domains.
%
%  Returns 1 value:
%
%    I - Scalar value of the integrated grid function.
%

% 13 May 2016
% Gustavo Rivera
% Cornell University

% Get some quadrature weights.
[junk, w, P] = lglnodes( double( n ) - 1 );

% Get some element widths.
Lx = max( x ) - min( x );
hx = Lx / mx;

for ii = 1:mx
    for jj = 1:n
        ndx    = (ii - 1) * n  + jj;
        u(ndx) = w(jj) .* u(ndx) .* (hx / 2);
    end
end

% Sum the kernel.
I = sum( u );

end
