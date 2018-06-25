function  [phi, psi] = flowfun( u, v, flag )
%  [phi, psi] = flowfun( u, v, flag )
%
% Computes the potential, phi, and the streamfunction, psi, of a 2-dimensional
% flow defined by the matrices of velocity components u and v, so that:
%
%           d(phi)    d(psi)          d(phi)    d(psi)
%      u =  -----  -  ----- ,    v =  -----  +  -----
%            dx        dy              dx        dy
%
% For a potential (irrotational) flow psi = 0, and the Laplacian of psi is
% equal to the divergence of the velocity field.  A non-divergent flow can be
% described by the streamfunction alone, and the Laplacian of the
% streamfunction is equal to vorticity (curl) of the velocity field.
% The stepsizes dx and dy are assumed to equal unity.
%
% Because these potentials are defined up to the integration constant their
% absolute values are such that phi(1, 1) = psi(1, 1) = 0.
%
% If only the streamfunction is needed, a flag may be supplied.  See below
% for how.
%
% Takes 3 arguments:
%
%  u    - 2D matrix of velocity components, either the difference of the
%         real-part of the derivatives, or the complex difference
%         (real-part) and summation (imaginary-part) of the derivatives.
%  v    - Optional 2D matrix of velocity components.  If u is specified
%         as only a difference of derivatives, v is the summation of
%         the derivatives.  Otherwise, not provided.
%  flag - Optional string to specify just the streamfunction, psi,
%         and not the potential, phi.  If specified as
%         ['-', 'psi', 'streamfunction'], then the output phi will
%         be omitted (returns psi only).
%
% Returns 2 values:
%
%  phi - 2D matrix, same size as u and v, containing the potential.
%        Will contain the streamfunction instead if flag was
%        supplied as such.
%  psi - 2D matrix, same size as u and v, containing the streamfunction.

%  Originally written by Kirill K. Pankratov, March 7, 1994.

% Check input arguments .............................................
issu   = 0;
issv   = 0;
isflag = 0;    % For input checking
isphi  = 1;
ispsi  = 1;        % Is phi and psi to be computed

if nargin == 1
    issu = isstr( u );
end
if nargin == 2
    issv = isstr( v );
end
if nargin == 1 & ~issu,
    v = imag( u );
end
if issv
    flag   = v;
    v      = imag( u );
    isflag = 1;
end
if nargin == 0 | issu            % Not enough input arguments
    disp( sprintf( '\n  Error: function must have input arguments:\n  matrivces  U and V  (or complex matrix W = U+iV)\n' ) )
  return
end
if any( size( u ) ~= size( v ) )     % Disparate sizes
    disp( sprintf( '\n  Error: matrices U and V must be of equal size\n' ) )
  return
end
if nargin == 3,
    isflag = 1;
end
u = real( u );

% Check the flag string . . . . . . . .
Dcn = char( '+', 'potential', 'phi' );
Dcn = char( Dcn, '-', 'streamfunction', 'psi' );
if isflag
  lmin = min( size( flag, 2 ), size( Dcn, 2 ) );
  flag = flag( 1, 1:lmin );
  A   = flag( ones( size( Dcn, 1 ), 1 ), 1:lmin ) == Dcn(:, 1:lmin);
  if lmin > 1
      coinc = sum( A' );
  else
      coinc = A';
  end
  fnd = find( coinc == lmin );
  if fnd ~= []
      if fnd < 4
          ispsi = 0;
      else
          isphi = 0;
      end end
end

phi = [];        % Create output
psi = [];

lx = size( u, 2 );  % Size of the velocity matrices
ly = size( u, 1 );

% Now the main computations .........................................
% Integrate velocity fields to get potential and streamfunction
% Use Simpson rule summation (function CUMSIMP)

 % Compute potential PHI (potential, non-rotating part)
if isphi
  cx  = cumsimp( u(1, :) );  % Compute x-integration constant
  cy  = cumsimp( v(:, 1) );  % Compute y-integration constant
  phi = cumsimp( v ) + cx(ones( ly, 1 ), :);
  phi = (phi + cumsimp( u' )' + cy(:, ones( 1, lx ))) / 2;
end

 % Compute streamfunction PSI (solenoidal part)
if ispsi
  cx  =  cumsimp( v(1, :) );  % Compute x-integration constant
  cy  =  cumsimp( u(:, 1) );  % Compute y-integration constant
  psi = -cumsimp( u ) + cx(ones( ly, 1 ), :);
  psi = (psi + cumsimp(v')' - cy(:, ones( 1, lx ))) / 2;
end

 % Rename output if need only PSI
 if ~isphi && ispsi && nargout==1
     phi = psi;
 end
