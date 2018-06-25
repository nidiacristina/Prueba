function [total_memory, memory_per_node] = smpm_check_mem_requirements( varargin )
% [total_memory, memory_per_node] = smpm_check_mem_requirements( n, mx, my, mz, np, ngmres_vis, ngmres_ppe, quiet_flag )
% [total_memory, memory_per_node] = smpm_check_mem_requirements( input_file_name, np, quiet_flag )
%
%
% Prints out the memory requirements for the SMPM based on the subdomains
% parameterization provided as arguments.  Requirements for total memory needed
% as well as per-node are reported in GiB.
%
% NOTE: This assumes that you're using deflated and preconditioned Schur, as
%       that is the only method that converges in reasonable time for 3D problems.
%
% NOTE: To maintain compatibility with previous versions, this function does not
%       return values unless explicitly requested.
%
% Takes 8 arguments:
%
%   n          - Number of GLL points per direction, per subdomain.
%   mx         - Number of subdomains in the x-direction.
%   my         - Number of transverse wavenumbers.
%   mz         - Number of subdomains in the z-direction.
%   np         - Number of processors used in the simulation.
%   ngmres_vis - Number of GMRES iterations for the viscous solvers.
%   ngmres_ppe - Number of pressure GMRES iterations for the SMW solver.
%   quiet_flag - Optional flag indicating no output should be generated.  If
%                omitted, defaults to false so output is generated.
%
% OR
%
% Takes 3 arguments:
%
%   input_file_name - Path of an input file.
%   np              - Number of processors used in the simulation.
%   quiet_flag      - Optional flag indicating no output should be generated.
%                     If omitted, defaults to false so output is generated.
%
% Returns 2 values:
%
%   total_memory    - Estimate of total number of bytes the specified problem
%                     requires across all MPI ranks.
%   memory_per_node - Estimate of number of bytes the specified problem requires
%                     for a single MPI rank.
%

% Parse the input arguments.
quiet_flag = false;

if nargin == 2 || nargin == 3
    data       = smpm_read_inputfile( varargin{1} );
    n          = data.n;
    mx         = data.nsubx;
    my         = data.nsuby;
    mz         = data.nsubz;
    np         = varargin{2};
    ngmres_vis = data.gmres_maxit_viscous;
    ngmres_ppe = data.gmres_maxit_poisson;

    if nargin == 3
        quiet_flag = varargin{3};
    end
elseif nargin >= 7
    n              = varargin{1};
    mx             = varargin{2};
    my             = varargin{3};
    mz             = varargin{4};
    np             = varargin{5};
    ngmres_vis     = varargin{6};
    ngmres_ppe     = varargin{7};
    if nargin > 7
        quiet_flag = varargin{8};
    end
end

% Establish a value for w.
w = 1;

% Set some constants.
r            = n * n * mx * mz;
subr         = n * n * mz * w;
subk         = 2 * n * mz * (w - 1);
subdimblock  = n * n * mz;
subnumblocks = w;
k            = 2 * n * mz * (mx - 1);
dimblock     = n * n * mz * w;
numblocks    = np;

% Estimate of the number of vectors of order r that we'll need.  Note that
% this should be an upper bound so that we don't actually count the number
% of long lived variables.
narrays = 40;

% Set some dimensions of the operator matrices.
s        = n * mz;
dimL     = [r, r];
dimC     = [k, k];
dimCk    = [4 * s, 2 * s, my];
dimMk    = [4 * s, 4 * s, my];
dimBk    = [k, dimblock];
dimsubAk = [subdimblock, subdimblock];
dimsubBk = [subr, subk];
dimsubC  = [subk, subk];

% Set some sizes for some of the storage arrays.
dimarrays = [r/np, narrays];

% Set some dimensions of the GMRES work arrays for the capacitance solver.
dimWorkV = [r/np, ngmres_vis];
dimWorkP = [k/np, ngmres_ppe];
       %
       % Note that both of the above assume that the
       % user is using parallel gmres, which uses
       % substantially less memory.

% Set some memory requirements in GiB.
doubles_per_GiB = 8 / 1024 / 1024 / 1024;
memL            = doubles_per_GiB * prod( dimL );
memC            = doubles_per_GiB * prod( dimC ) * my * 2;
memCk           = doubles_per_GiB * prod( dimCk );
memMk           = doubles_per_GiB * prod( dimMk );
memBk           = 0 * doubles_per_GiB * prod( dimBk );
memsubC         = 0 * doubles_per_GiB * prod( dimsubC );
memsubAk        = doubles_per_GiB * prod( dimsubAk ) * 2; % The Schur factorizations require more storage.
memsubBk        = 0 * doubles_per_GiB * prod( dimsubBk );
memarrays       = doubles_per_GiB * prod( dimarrays ) * my;
memgmres        = doubles_per_GiB * (prod( dimWorkV ) + prod( dimWorkP )) * my;

% Compute the memory requirements per node and in total.
memory_per_node  = ((mx / np) * memsubAk + ...
                    (mx / np) * memMk + ...
                    (mx / np) * memCk + ...
                    memgmres + memarrays);
total_memory  = memory_per_node * np;

% Print our requirements if the user has requested them.
if ~quiet_flag
    fprintf( '\n');
    fprintf( '  Simulation total memory (GiB) : %8.3f \n', total_memory );
    fprintf( '\n');
    fprintf( '                   n : %8d \n', n );
    fprintf( '               nsubx : %8d \n', mx );
    fprintf( '               nsuby : %8d \n', my );
    fprintf( '               nsubz : %8d \n', mz );
    fprintf( '  degrees of freedom : %8d \n', r * my );
    fprintf( '               ranks : %8d \n', np );
    fprintf( '         Ak per rank : %8d \n', mx / np );
    fprintf( '             dim(Ak) : %8d \n', dimblock );
    fprintf( '          dim(subAk) : %8d \n', dimsubAk(1) );
    fprintf( '             dim(Sk) : %8d, %8d, %8d \n', dimCk(1), dimCk(2), dimCk(3) );
    fprintf( '             dim(Mk) : %8d, %8d, %8d \n', dimMk(1), dimMk(2), dimMk(3) );

    fprintf( ' \n');
    fprintf( '  Memory per process (GiB) : %8.3f \n', memory_per_node );
    fprintf( ' \n');
    fprintf( '     local Ak blocks : %8.3f \n', (mx / np) * memsubAk );
    fprintf( '     local Sk blocks : %8.3f \n', (mx / np) * memCk );
    fprintf( '     local Mk blocks : %8.3f \n', (mx / np) * memMk );
    fprintf( '         data arrays : %8.3f \n', memarrays );
    fprintf( '   gmres work arrays : %8.3f \n', memgmres );
end

% Convert the requirements into bytes for the caller so they may be used
% programmatically.
if nargout > 0
    total_memory    = total_memory * 1024^3;
    memory_per_node = memory_per_node * 1024^3;
else
    clear total_memory memory_per_node;
end

end
