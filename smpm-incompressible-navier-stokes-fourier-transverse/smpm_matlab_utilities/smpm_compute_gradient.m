function [ phi_x, phi_y, phi_z ] = smpm_compute_gradient( phi, n, mx, my, mz, x, y, z )
% [ phi_x, phi_z  ] = smpm_compute_gradient( phi, n, mx, mz, x, z )
%
% Computes the gradient of a function, phi, defined in a 2D GLL grid with a
% transverse, equally spaced, domain. The differentiation is performed in
% 2D via Spectral Differentiation Matrix based on GLL quadrature and in the
% transverse via the DFT.
%
% Takes 8 arguments:
%
%   n   - Number of GLL points per direction, per subdomain.
%   mx  - Number of subdomains in the x-direction.
%   my  - Number of points in the y-direction.
%   mz  - Number of subdomains in the z-direction.
%   x   - Column vector, of size [n * mz, n * mx , my], containing the x
%         coordinates of the generated mesh.
%   y   - Column vector, of size [n * mz, n * mx , my], containing the y
%         coordinates of the generated mesh.
%   z   - Column vector, of size [n * mz, n * mx , my], containing the z
%         coordinates of the generated mesh.
%   phi - Column vector, of size [n * mz, n * mx , my], containing the
%         function to be differentiated. 
%
% Returns 3 values:
%
%   phi_x - derivative of the function, of size [n * mz, n * mx , my], in 
%           the x-direction.
%   phi_y - derivative of the function, of size [n * mz, n * mx , my], in 
%           the y-direction.
%   phi_z - derivative of the function, of size [n * mz, n * mx , my], in 
%           the z-direction.
%
% 16 August 2017
% Gustavo Rivera

%% Prepare Arrays and compute some relevant parameters

% Reshape the field into a 2D vector
phi_in = reshape(phi, n * n * mx * mz, my);

% Initalize 2D derivative arrays
vx = zeros( n * mz * n * mx, my);
vz = zeros( n * mz * n * mx, my);
vy = zeros( n * mz * n * mx, my);

% Initialize 3D derivative arrays
phi_x =  zeros( n * mz, n * mx, my);
phi_y =  zeros( n * mz, n * mx, my);
phi_z =  zeros( n * mz, n * mx, my);


%% Step 1: Build the Deformation Maps

% Compute GLL Nodes
xi = lglnodes( double( n ) - 1 );

% Compute the differentiation matrix
[D, ~, ~] = derv( n - 1, xi , n - 1 );
D         = reshape( D, n, n);

% Reshape arrays into 2D
xin = reshape(x, n * mz * n * mx, my );
zin = reshape(z, n * mz * n * mx, my );

% Set the field by taking the first transverse point
iix = xin(:,1);
iiz = zin(:,1);

% Compute derivative in xi
iix_xi = D * reshape( iix, n , n * mx * mz );
iiz_xi = D * reshape( iiz, n , n * mx * mz );

% Permute into x/eta-first indexing.
iix    = perfect_shuffle( n*mx, n*mz, iix);
iiz    = perfect_shuffle( n*mx, n*mz, iiz);

iix_xi = iix_xi(:); 
iiz_xi = iiz_xi(:); 

% Compute derivative in eta
iix_eta = D * reshape( iix, n , n * mx * mz );
iiz_eta = D * reshape( iiz, n , n * mx * mz );

% Permute back into xi/z-indexing.
iix     = perfect_shuffle( n*mz, n*mx, iix);
iiz     = perfect_shuffle( n*mz, n*mx, iiz);

iix_eta = iix_eta(:); 
iiz_eta = iiz_eta(:); 

% Build the determinant of the Jacobian matrix and its derivatives.
iidetJ  = iix_eta .* iiz_xi - iiz_eta .* iix_xi;

% Store the data in the global vectors.
x_xi   = iix_xi;
z_xi   = iiz_xi;
x_eta  = iix_eta;
z_eta  = iiz_eta;
detJ   = iidetJ;

%% Step 2: Differentiate field in the x-z plane

% Loop per transverse plane
for ii = 1:my

    % Store the transverse plane
    iiv = phi_in(:,ii);

    % Compute Derivative in xi
    iiv_xi = D * reshape(iiv,n,n*mx*mz);
    iiv_xi = iiv_xi(:);

    % Shuffle density into eta-first index
    iiv    = perfect_shuffle( n*mx, n*mz, iiv);

    % Compute Derivative in eta
    iiv_eta =  D * reshape(iiv,n,n*mx*mz);
    iiv_eta = iiv_eta(:);

    % Shuffle back into xi-first indexing
    iiv_eta = perfect_shuffle( n*mz, n*mx, iiv_eta);

    % Map into physical domain
    vx(:,ii) = ( z_xi .* iiv_eta' - z_eta .* iiv_xi )./ detJ;
    vz(:,ii) = ( x_eta .* iiv_xi - x_xi .* iiv_eta' )./ detJ;
end

% Reshape Derivatives in x-z
phi_x = reshape( vx, n * mz, n * mx, my );
phi_z = reshape( vz, n * mz, n * mx, my );

%% Step 3: Differentiate field in the y plane


% We will only differentiate if there is a transverse domain
if my > 1

    % Compute the grid spacing in the transverse and determine transverse
    % domain length
    dy = y(1,1,2) - y(1,1,1);
    Ly = max(y(:)) - min(y(:)) + dy;

    % Build wavenumbers
    wavenumber = 2*pi/(Ly).*[0:my/2-1, 0, -my/2+1:-1];

    % Global number of GLL Points
    nsg = n * n * mx * mz;

    % FFT each row of the field
    for ii = 1:nsg
        phif(ii,:) = fft( phi_in(ii,:) );
    end
    
    % Compute fourier derivative
    for jj = 1:my
        dphif(:,jj) = 1i .* wavenumber(jj) .* phif(:,jj);
    end
    
    % IFFT each row of the derivative
    for ii = 1:nsg
        vy(ii,:) = ifft( dphif(ii,:) );
    end
    
    % Reshape Derivative in y
    phi_y = reshape( vy, n * mz, n * mx, my );

end
end

