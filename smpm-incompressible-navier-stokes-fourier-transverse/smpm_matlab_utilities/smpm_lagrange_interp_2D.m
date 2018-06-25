function qi = smpm_lagrange_interp_2D( x, z, q, xi, zi );
% qi = smpm_lagrange_interp_2D( x, z, q, xi, zi );
%
% Computes the single element 2D Lagrange interpolant of q given
% GLL points x and z, and interpolation points xi, zi.
%
% Takes 5 arguments:
%
%   x  - Vector, length (n x 1), of GLL point x-coordinates.
%   z  - Vector, length (n x 1), of GLL point z-coordinates.
%   q  - Matrix, sized (n x n), of the grid function to interpolate.
%   xi - Vector of interpolation point x-coordinates.
%   zi - Vector of interpolation point z-coordinates.
%
% Returns 1 value:
%
%   qi  - The interpolant.
%

% Make sure all inputs are column vectors.
x  = x(:);
z  = z(:);
xi = xi(:);
zi = zi(:);

% Get some constants.
n = length( x );

% Construct the 1D Lagrange interpolant in x.
Lx = @(w,k)(prod( (w - x(x ~= x(k))) ./(x(k) - x(x ~= x(k))) ));

% Construct the 1D Lagrange interpolant in z.
Lz = @(w,k)(prod( (w - z(z ~= z(k))) ./ (z(k) - z(z ~= z(k))) ));

% Compute all the Lagrange polynomials at all the interpolation points.
Lxi = zeros( length(xi), n );
for ii = 1:n
    for jj = 1:length( xi )
        Lxi(jj, ii) = Lx( xi(jj), ii );
        Lzi(jj, ii) = Lz( zi(jj), ii );
    end
end

% Compute the interpolation.
qi = 0 * xi;
for ii = 1:n
    for jj = 1:n
        qi =  qi + q(jj, ii) * Lxi(:, ii) .* Lzi(:, jj);
    end
end

end
