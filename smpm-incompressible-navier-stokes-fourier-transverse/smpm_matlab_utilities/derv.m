function [d, d2, d3] = derv( nterm, x, ndim )
% [d, d2, d3] = derv( nterm, x, ndim )
%
% Calculate the first, second, and third derivative matrix at
% Gauss-Lobatto-Legendre grid points calculated in the Fortran routines
% jacobl() and jacobf().
%
% The theory behind this code can be found in Costa & Don "On the computation
% of high order derivatives".  Applied Numerical Mathematics. 33 (2000) 151 -
% 159.
%
% Takes 3 arguments:
%
%   nterm - Polynomial degree.
%   x     - Gauss-Lobatto-Legendre Points.
%   ndim  - Number of grid points in each direction.
%
% Returns 3 values:
%
%   d  - Vector, length (n+1)*(n+1), of first derivatives in column-wise order.
%   d2 - Vector, length (n+1)*(n+1), of second derivatives in column-wise order.
%   d3 - Vector, length (n+1)*(n+1), of third derivatives in column-wise order.

% Adapted from Pete Diamessis' derv.f and Sumedh Joshi's derv.f90.
%
% Gustavo Rivera 2016

% Allocate derivative matrices.
d  = zeros( ndim+1, ndim+1 );
d2 = zeros( ndim+1, ndim+1 );
d3 = zeros( ndim+1, ndim+1 );

% 1) First compute coefficients C_i.  See Eqn (4) in Costa and Don.
for k = 1: nterm+1
    prod = 1.;
    xk   = x(k);

    for l = 1:nterm+1
        xl = x(l);

        if l ~= k
            prod = prod * (xk - xl);
        else
            prod = prod;
        end
    end

    c(k) = prod;
end

%% Compute the first derivative.

% Calculate off-diagonal elements of D^1.  See Eqn (6).
for k = 1: nterm+1
    xk = x(k);

    for j = 1:nterm+1
        xj = x(j);

        if k ~= j
            d(k, j) = c(k) / (c(j) * (xk - xj));
        else
            d(k, j) = 0.;
        end
    end
end

% Calculate diagonal elements now.  Diagonal element is negative row sum of
% off diagonal elements.  See Costa & Don Eqn(9).
for k = 1: nterm+1
    sum = 0.;

    for j = 1: nterm+1
        if k ~= j
            sum = sum + d(k, j);
        else
            sum = sum;
        end
    end

    d(k, k) = -sum;
end

%% Compute the second derivative.

% Off-diagonal elements.  See Eqn(13).
for k = 1: nterm+1
    xk = x(k);

    for j = 1: nterm+1
        xj = x(j);

        if k ~= j
            d2(k, j) = 2. * (d(k, k) * d(k, j) - d(k, j) / (xk - xj));
        else
            d2(k, j) = 0.;
        end
    end
end

% Diagonal elements.  See Eqn(9).
for k = 1: nterm+1
    sum = 0.;

    for j = 1:nterm+1

        if k ~= j
            sum = sum + d2(k, j);
        else
            sum = sum;
        end
    end

    d2(k, k) = -sum;
end

%% Compute the third derivative.

% Off-diagonal elements.
for k = 1: nterm+1
    xk = x(k);

    for j = 1: nterm+1
        xj = x(j);

        if k ~= j
            d3(k, j) = 3. * (d2(k, k) * d(k, j) - d2(k, j) / (xk - xj));
        else
            d3(k, j) = 0.;
        end
    end
end

% Diagonal elements.
for k = 1: nterm+1
    sum = 0.;

    for j = 1: nterm+1
        if k ~= j
            sum = sum + d3(k, j);
        else
            sum = sum;
        end
    end

    d3(k, k) = -sum;
end

% Reverse the ordering of the matrices and flatten them into vectors.
d  = flipud(  d(:) );
d2 = flipud( d2(:) );
d3 = flipud( d3(:) );

end
