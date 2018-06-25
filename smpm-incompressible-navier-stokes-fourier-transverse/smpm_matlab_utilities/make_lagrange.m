function y = make_lagrange( x, xn, n, d )
% y = make_lagrange( x, xn, n, d )
%
% Computes the Lagrange polynomial l_n(x) or its d-th derivative.
%
% Takes 4 arguments:
%
%   x  - Row vector of x-coordinates to evaluate.
%   xn - Row vector of knots.
%   n  - Index, in the range of [1, length( x )], of l_n(x).
%   d  - Optional derivative order to compute.  If omitted, defaults to 0.
%
% Returns 1 value:
%
%   y  - The computed Lagrange polynomial, l_n(x).
%

% Sumedh Joshi
% 9/20/2011

if nargin < 4
    d = 0;
end

% Compute the coefficients of the polynomial in the canonical form.
x0 = xn(n);
xn = [xn(1:n-1), xn(n+1:end)];
N  = length( xn );
An = zeros( 1, N+1 ); % An = [a(0) a(1) a(2) ... a(N)];

for ii = 0:N
    % Get the set of all subsets of xn of length ii+1.
    iiT = nchoosek( xn, ii );

    % Compute the sum of the products of each subset.
    An(ii+1) = (-1)^(ii)*sum( prod( iiT, 2 ) );
end

% Compute the denominator.
denominator = prod( x0 - xn );

% Compute the coefficients.
An = An / denominator;

% Compute the derivatives.
if d == 0
    % Compute the polynomial.
    y = polyval( An, x );
else
    for ii = 1:d
        iiN = length( An );
        c   = [iiN:-1:1];
        c   = c - 1;
        An  = An .* c;
        An  = An(1:end-1);
    end
    y = polyval( An, x );
end

end
