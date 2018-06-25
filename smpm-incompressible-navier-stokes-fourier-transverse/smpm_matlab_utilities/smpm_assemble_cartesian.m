function [A, E, B, S] = smpm_assemble_cartesian( varargin )
% [A, E, B, S] = smpm_assemble_cartesian( n, mx, mz, Lx, Lz[, t] )
%
%  Builds the domain decomposition of the spectral multidomain penalty method
%  discretization of the Poisson equation on a uniform cartesian grid.
%  Returns all matrices as sparse objects.
%
%  To assemble the SMPM matrix, compute:
%
%     SMPM = A + E * B;
%
%  Takes 5 or 6 arguments:
%
%    n  - Number of GLL points per direction, per subdomain.
%    mx - Number of subdomains in the x-direction.
%    mz - Number of subdomains in the z-direction.
%    Lx - Length between the first and last elements in the x-direction.
%    Lz - Length between the first and last elements in the z-direction.
%    t  - Optional parameter indicating how many vertical strips are in
%         one sub-region in the domain decomposition.  Must evenly divide
%         mx and defaults to 1 if omitted.
%
%  Returns 4 values:
%
%    A -
%    E -
%    B -
%    S -
%
%  11 Feb 2015
%  Sumedh Joshi

% Parse the arguments.
if nargin == 5
   t = 1;
   n  = varargin{1};
   mx = varargin{2};
   mz = varargin{3};
   Lx = varargin{4};
   Lz = varargin{5};
else
   n  = varargin{1};
   mx = varargin{2};
   mz = varargin{3};
   Lx = varargin{4};
   Lz = varargin{5};
   t  = varargin{6};
end

% Set the polynomial order in each element in each direction.
p = n - 1;

% Set the patching penalty factor.
F = 1000;

% Compute the element-wise coordinates.

% Compute the GLL gridpoints.
n            = p;
[x, w, junk] = lglnodes( p );
x            = flipud( x )';

% Compute the spectral differentiation matrix of order 2.
D2 = zeros( n+1, n+1 );
for ii = 1:n+1
    D2(:, ii) = make_lagrange( x, x, ii, 2 );
end

% Compute the spectral differentiation matrix of order 1.
D1 = zeros( n+1, n+1 );
for ii = 1:n+1
    D1(:, ii) = make_lagrange( x, x, ii, 1 );
end
x = x';

% Build the grid in xi-first indexing.
n      = n + 1;
r      = n * n * mx * mz;
[x, z] = smpm_build_cartesian_mesh( n, mx, mz, [0, Lx], [0, Lz] );

% Set some deformation constants (in place of maps since we assume the grid is plaid).
hx = Lx / mx;
hz = Lz / mz;

%% Build a 1D SMPM matrix in x.

% The 1st derivative in x.
Dx  = (2 / hx) * kron( speye( mx, mx ), D1 );

% The 2nd derivative in x.
Dxx = (2 / hx)^2 * kron( speye( mx, mx ), D2 );

% The C0 continuity conditions.
C0x = sparse( Dxx * 0 );
for ii = 1:mx - 1
    left                        = ii * n;
    right                       = ii * n + 1;
    C0x(left:right, left:right) = [1 -1; -1 1];
end

% The C1 continuity conditions.
C1x = sparse( Dxx * 0 );
for ii = 1:mx - 1
    left                        = ii * n;
    right                       = ii * n + 1;
    C1x(left:right, left:right) = [ 1 -1; 1 -1];
end

% The Neumann boundary conditions.
B1x           = sparse( Dxx * 0 );
B1x(1,   1)   = -1;
B1x(end, end) =  1;

% Set the penalty parameter.
w0     = w(1);
taumin = (1 / w0) * (1 + 2*w0 - 2*sqrt( w0^2 + w0 )) * 2 / hx;
taumax = (1 / w0) * (1 + 2*w0 + 2*sqrt( w0^2 + w0 )) * 2 / hx;
tau    = (taumax + taumin) * 0.25;

% Assemble the matrix.
SMPMx = Dxx + tau * F * (C0x + sparse( C1x * Dx )) + tau * sparse( B1x * Dx );

%% Build a 1D SMPM matrix in z.

% The 1st derivative in z.
Dz  = (2 / hz) * kron( speye( mz, mz ), D1 );

% The 2nd derivative in x.
Dzz = (2 / hz)^2 * kron( speye( mz, mz ), D2 );

% The C0 continuity conditions.
C0z = sparse( Dzz * 0 );
for ii = 1:mz - 1
    left                        = ii * n;
    right                       = ii * n + 1;
    C0z(left:right, left:right) = [1 -1; -1 1];
end

% The C1 continuity conditions.
C1z = sparse( Dzz * 0 );
for ii = 1:mz - 1
    left                        = ii * n;
    right                       = ii * n + 1;
    C1z(left:right, left:right) = [ 1 -1; 1 -1];
end

% The Neumann boundary conditions.
B1z           = sparse( Dzz * 0 );
B1z(1,   1)   = -1;
B1z(end, end) =  1;

% Set the penalty parameter.
w0     = w(1);
taumin = (1 / w0) * (1 + 2*w0 - 2*sqrt( w0^2 + w0 )) * 2 / hz;
taumax = (1 / w0) * (1 + 2*w0 + 2*sqrt( w0^2 + w0 )) * 2 / hz;
tau    = (taumax + taumin) * 0.25;

% Assemble the matrix.
SMPMz = Dzz + tau * F * (C0z + sparse( C1z * Dz )) + tau * sparse( B1z * Dz );

%% Build the 2D SMPM matrix from the 1D SMPM matrix.
SMPM = kron( speye( mx * n ), ...
             sparse( SMPMz ) ) + ...
       kron( sparse( SMPMx ), ...
             speye( mz * n ) );

%% Assemble the Schur matrix and the other components of the domain decomposition.

% Set some constants we'll need.
s = mz * n;
k = 2 * s * (mx / t - 1);

% Grab the block-diagonal internal blocks.
A = SMPM .* kron( speye( mx / t , mx / t ), ...
                  ones( n * n * mz * t, n * n * mz * t ) );

% Grab the off-diagaonl pre-Schur matrix.
EB = sparse( SMPM - A );

% Build the extension operator.
E   = sparse( r, k );
ndx = sum( abs( EB ), 2 );
ndx = find( ndx ~= 0 );
for ii = 1:length( ndx )
    E(ndx(ii), ii) = 1;
end

% Build the off-block flux matrix B.
B = E' * EB;

% Build the Schur matrix.
dimblock = n * n * mz * t;
S        = speye( k );
for ii = 1:mx
    a = (ii - 1) * dimblock + 1;
    z = (ii - 1) * dimblock + dimblock;
    S = S + B(:, a:z) * sparse( A(a:z, a:z) \ E(a:z, :) );
end

end
