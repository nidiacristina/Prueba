function smpm_build_test_cases( test_root_directory )
% smpm_build_test_cases( test_root_directory )
%
% Creates a number of test cases for the SMPM incompressible Navier-Stokes
% solver.  This creates input configurations for test cases that exercise
% different aspects of the code base suitable for correctness and performance
% testing.
%
% Takes 1 argument:
%
%   testing_root_directory - Optional string specifying the directory to
%                            generate the test cases in.  If omitted, defaults
%                            to the current directory.
%
% Returns nothing.

if nargin < 1
   test_root_directory = [];
end

if isempty( test_root_directory )
   test_root_directory = pwd;
end

% Enumerate the function names that will build our test cases.  These will
% be created in the order specified in directories beneath test_root_directory.
test_case_builders = { 'build_smpm_test_lid_driven_cavity_undeformed', ...
                       'build_smpm_test_lid_driven_cavity_deformed', ...
                       'build_smpm_test_boussinesq_lock_exchange', ...
                       'build_smpm_test_boussinesq_two_layer' };

% Build each of the test cases.
for test_case_number = 1:numel( test_case_builders )
    test_directory = sprintf( '%s/%d', ...
                              test_root_directory, test_case_number );
    mkdir( test_directory );

    feval( test_case_builders{test_case_number}, test_directory );
end

return

function [configuration_file_name, conditions_file_name, ...
          tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name )
% [configuration_file_name, conditions_file_name, ...
%          tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name )
%
% Builds file names from the supplied test directory and name.
%
% Takes 2 arguments:
%
%   test_directory - Directory where the test is to be built.
%   test_name      - Name of the test being built.
%
% Returns 4 values:
%
%   configuration_file_name - Path to the test's configuration file.
%   conditions_file_name    - Path to the test's initial conditions file.
%   tolerance_file_name     - Path to the test's tolerance file.
%   gold_file_name          - Path to the test's gold output file.
%

   configuration_file_name = sprintf( '%s/%s_in', ...
                                      test_directory, test_name );
   conditions_file_name    = sprintf( '%s/%s_init.h5', ...
                                      test_directory, test_name );
   tolerance_file_name     = sprintf( '%s/%s_tolerance', ...
                                      test_directory, test_name );
   gold_file_name          = sprintf( '%s/%s_gold', ...
                                      test_directory, test_name );

return

function build_smpm_test_lid_driven_cavity_undeformed( test_directory )
% build_smpm_test_lid_driven_cavity_undeformed( test_directory )
%
% XXX:
%
% Takes 1 argument:
%
%   test_directory - Directory where the test is to be built.
%
% Returns nothing.

   % Test Case: lid-driven cavity (undeformed).
   test_name = 'lid_driven_cavity';

   [configuration_file_name, conditions_file_name, ...
    tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name );

   % Delete old versions of the run files.
   delete(  configuration_file_name );
   delete(  conditions_file_name );

   % Set some constants.
   n     = 8;
   nsubx = 16;
   nsubz = 16;
   Lx    = 1.0;
   Lz    = 1.0;

   % Build the mesh.
   [x, z] = smpm_build_cartesian_mesh( n, nsubx, nsubz, [0 Lx], [0 Lz] );

   % Set the initial condition for the lid-driven cavity.
   ux0          = 0*x;
   uz0          = 0*x;
   rho0         = 0*x;
   rho0_bar     = 0*x;
   s            = n * nsubz;
   ux0(s:s:end) = 16 * x(s:s:end).^2 .* (1 - x(s:s:end).^2);

   % Write the initial conditions file.
   smpm_write_initfile( n, nsubx, nsubz, x, z, rho0, ux0, uz0, rho0_bar, conditions_file_name );

   % Build the inputs for the input file.
   inputs = build_smpm_test_inputs( test_name );

   inputs.n                              = n;
   inputs.nsubx                          = nsubx;
   inputs.nsubz                          = nsubz;
   inputs.tend                           = 0.03;
   inputs.dt                             = 3e-4;
   inputs.facrobin                       = 10000;
   inputs.facrobin_ppe                   = 100;
   inputs.nu                             = 0.0040;
   inputs.nu_d                           = 0.0040;
   inputs.filter_order_xz                = 11;
   inputs.filter_order_y                 = 11;
   inputs.gmres_tol_poisson              = 1e-6;
   inputs.gmres_maxit_poisson            = 1000;
   inputs.gmres_tol_viscous              = 1e-6;
   inputs.gmres_maxit_viscous            = 200;
   inputs.bc_viscous_x                   = [1, 1, 1, 1];
   inputs.bc_viscous_z                   = [1, 1, 1, 1];
   inputs.bc_diffusion                   = [1, 1, 1, 1];
   inputs.check_null_error               = logical(1);
   inputs.check_numerical_error          = logical(1);
   inputs.read_bcs_from_initfile         = logical(1);
   inputs.use_capacitance_preconditioner = logical(1);
   inputs.use_deflation                  = logical(1);

   % Write the input file.
   smpm_write_inputfile( configuration_file_name, inputs );

   % Write the tolerance file.
   fid = fopen( tolerance_file_name, 'w+' );
   fprintf( fid, '1e-10' );
   fclose( fid );

   % Write a dummy gold out file.
   fid = fopen( gold_file_name, 'w+' );
   fprintf( fid, 'This file is intentionally empty.' );
   fclose( fid );

return

function build_smpm_test_lid_driven_cavity_deformed( test_directory )
% build_smpm_test_lid_driven_cavity_deformed( test_directory )
%
% XXX:
%
% Takes 1 argument:
%
%   test_directory - Directory where the test is to be built.
%
% Returns nothing.

   % Test Case: lid-driven cavity (internally deformed interfaces).
   test_name = 'lid_driven_cavity';

   [configuration_file_name, conditions_file_name, ...
    tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name );

   % Delete old versions of the run files.
   delete(  configuration_file_name );
   delete(  conditions_file_name );

   % Set some constants.
   n     = 15;
   nsubx = 10;
   nsubz = 10;
   Lx    = 1.0;
   Lz    = 1.0;

   % Construct the mesh.
   nint           = 1024;
   xlayer         = linspace( 0, Lx, nint );
   zlayer         = zeros( nsubz + 1, nint );
   zlayer(1, :)   = 0;
   zlayer(end, :) = 1;

   for ii = 1:nsubz - 1
       % Get the mean of this z-layer.
       iiz0            = 0 + ii * Lz / nsubz;

       % Build the layer with a sinusoidal oscillation.
       zlayer(ii+1, :) = iiz0 + sign( (-1)^ii ) * sin( xlayer * 10 ) / 100;
       %zlayer(ii+1, :) = iiz0 + sign( (-1)^ii ) * ( xlayer/50 - ( xlayer - 0.5).^3/10 + xlayer.^2/100 )
   end

   % Build the mesh.
   [x, z] = smpm_build_layered_mesh( n, nsubx, nsubz, xlayer, zlayer );

   % Perfect shuffle the mesh so that we have a sinusoidally deformed x-coordinate too.
   x2      = 0 * x;
   s       = n * nsubz;
   counter = 1;

   for jj = 1:s
       for kk = 1:s
           ndx         = (kk-1) * s + jj;
           x2(counter) = z(ndx);
           counter     = counter + 1;
       end
   end
   x = x2;

   % Set the initial condition for the lid-driven cavity.
   ux0          = 0*x;
   uz0          = 0*x;
   rho0         = 0*x;
   rho0_bar     = 0*x;
   s            = n * nsubz;
   ux0(s:s:end) = 16 * x(s:s:end).^2 .* (1 - x(s:s:end).^2);

   % Write the initial conditions file.
   smpm_write_initfile( n, nsubx, nsubz, x, z, rho0, ux0, uz0, rho0_bar, conditions_file_name );

   % Build the inputs for the input file.
   inputs = build_smpm_test_inputs( test_name );

   inputs.n                              = n;
   inputs.nsubx                          = nsubx;
   inputs.nsubz                          = nsubz;
   inputs.tend                           = 0.03;
   inputs.dt                             = 3e-4;
   inputs.facrobin                       = 10000;
   inputs.facrobin_ppe                   = 100;
   inputs.nu                             = 0.0040;
   inputs.nu_d                           = 0.0040;
   inputs.filter_order_xz                = 11;
   inputs.filter_order_y                 = 11;
   inputs.gmres_tol_poisson              = 1e-6;
   inputs.gmres_maxit_poisson            = 1000;
   inputs.gmres_tol_viscous              = 1e-6;
   inputs.gmres_maxit_viscous            = 200;
   inputs.bc_viscous_x                   = [1, 1, 1, 1];
   inputs.bc_viscous_z                   = [1, 1, 1, 1];
   inputs.bc_diffusion                   = [1, 1, 1, 1];
   inputs.check_null_error               = logical(1);
   inputs.check_numerical_error          = logical(1);
   inputs.read_bcs_from_initfile         = logical(1);
   inputs.use_capacitance_preconditioner = logical(1);
   inputs.use_deflation                  = logical(1);

   % Write the input file.
   smpm_write_inputfile( configuration_file_name, inputs );

   % Write the tolerance file.
   fid = fopen( tolerance_file_name, 'w+' );
   fprintf( fid, '1e-10' );
   fclose( fid );

   % Write a dummy gold out file.
   fid = fopen( gold_file_name, 'w+' );
   fprintf( fid, 'This file is intentionally empty.' );
   fclose( fid );

return

function build_smpm_test_boussinesq_lock_exchange( test_directory )
% build_smpm_test_boussinesq_lock_exchange( test_directory )
%
% XXX:
%
% Takes 1 argument:
%
%   test_directory - Directory where the test is to be built.
%
% Returns nothing.

   % Test Case: Boussinesq stratified lock exchange problem.
   test_name = 'lock_exchange';

   [configuration_file_name, conditions_file_name, ...
    tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name );

   % Delete old versions of the run files.
   delete(  configuration_file_name );
   delete(  conditions_file_name );

   % Set some constants (these are Jorge's values).
   n       = 15;
   nsubx   = 60;
   nsubz   = 10;
   xlim    = [0, 0.8];
   zlim    = [0, 0.1];

   % Build the mesh file.
   Z      = linspace( zlim(1), zlim(2), nsubz + 1 );
   X      = linspace( -1, 1, nsubx + 1 );
   X      = (X - X(1)) / X(end) * xlim(2) / 2;
   [x, z] = smpm_build_kron_mesh( n, X, Z );

   % Build the initial conditions file.
   espi       = 1000;
   rho        = 0.5 * (1 + tanh( espi * (x - 0.4) )) - 0.5;
   ux         = 0;
   uz         = 0;
   rho_bar    = 0*rho;
   smpm_write_initfile( n, nsubx, nsubz, x, z, rho, ux, uz, rho_bar, conditions_file_name );

   % Build the inputs for the input file.
   inputs = build_smpm_test_inputs( test_name );

   inputs.n                              = n;
   inputs.nsubx                          = nsubx;
   inputs.nsubz                          = nsubz;
   inputs.tend                           = 1.0;
   inputs.dt                             = 1e-3;
   inputs.facrobin                       = 10000;
   inputs.facrobin_ppe                   = 10000;
   inputs.nu                             = 1e-6;
   inputs.nu_d                           = 1e-6;
   inputs.filter_order_xz                = 11;
   inputs.filter_order_y                 = 11;
   inputs.gmres_tol_poisson              = 1e-8;
   inputs.gmres_maxit_poisson            = 1000;
   inputs.gmres_tol_viscous              = 1e-8;
   inputs.gmres_maxit_viscous            = 100;
   inputs.bc_viscous_x                   = [2, 1, 2, 1];
   inputs.bc_viscous_z                   = [1, 2, 1, 2];
   inputs.bc_diffusion                   = [2, 2, 2, 2];
   inputs.check_null_error               = logical(0);
   inputs.check_numerical_error          = logical(1);
   inputs.read_bcs_from_initfile         = logical(0);
   inputs.use_capacitance_preconditioner = logical(1);
   inputs.use_deflation                  = logical(1);
   inputs.timesteps_between_writes       = 100;
   inputs.do_interfacial_averaging       = logical(0);

   % Write the input file.
   smpm_write_inputfile( configuration_file_name, inputs );

   % Write the tolerance file.
   fid = fopen( tolerance_file_name, 'w+' );
   fprintf( fid, '1e-10' );
   fclose( fid );

   % Write a dummy gold out file.
   fid = fopen( gold_file_name, 'w+' );
   fprintf( fid, 'This file is intentionally empty.' );
   fclose( fid );

return

function build_smpm_test_boussinesq_two_layer( test_directory )
% build_smpm_test_boussinesq_two_layer( test_directory )
%
% XXX:
%
% Takes 1 argument:
%
%   test_directory - Directory where the test is to be built.
%
% Returns nothing.

   % Test Case: Boussinesq stably-stratified 2-layer fluid.
   test_name = 'two_layer';

   [configuration_file_name, conditions_file_name, ...
    tolerance_file_name, gold_file_name] = get_file_names( test_directory, test_name );

   % Delete old versions of the run files.
   delete(  configuration_file_name );
   delete(  conditions_file_name );

   % Set some constants (these are Jorge's values).
   n       = 5;
   nsubx   = 20;
   nsubz   = 10;
   xlim    = [0, 0.8];
   zlim    = [0, 0.1];

   % Build the mesh file.
   Z      = linspace( zlim(1), zlim(2), nsubz + 1 );
   X      = linspace( -1, 1, nsubx + 1 );
   X      = (X - X(1)) / X(end) * xlim(2) / 2;
   [x, z] = smpm_build_kron_mesh( n, X, Z );

   % Build the initial conditions file.
   espi       = 500;
   rho        = -1.0 * (0.5 * (1 + tanh( espi * (z - 0.05) )) - 0.5);
   ux         = 0;
   uz         = 0;
   rho_bar    = 0*rho;
   smpm_write_initfile( n, nsubx, nsubz, x, z, rho, ux, uz, rho_bar, conditions_file_name );

   % Build the inputs for the input file.
   inputs = build_smpm_test_inputs( test_name );

   inputs.n                              = n;
   inputs.nsubx                          = nsubx;
   inputs.nsubz                          = nsubz;
   inputs.tend                           = 1.0;
   inputs.dt                             = 0.01;
   inputs.facrobin                       = 10000;
   inputs.facrobin_ppe                   = 10000;
   inputs.nu                             = 1e-6;
   inputs.nu_d                           = 1e-6;
   inputs.filter_order_xz                = 11;
   inputs.filter_order_y                 = 11;
   inputs.gmres_tol_poisson              = 1e-8;
   inputs.gmres_maxit_poisson            = 1000;
   inputs.gmres_tol_viscous              = 1e-8;
   inputs.gmres_maxit_viscous            = 100;
   inputs.bc_viscous_x                   = [2, 1, 2, 1];
   inputs.bc_viscous_z                   = [1, 2, 1, 2];
   inputs.bc_diffusion                   = [2, 2, 2, 2];
   inputs.check_null_error               = logical(0);
   inputs.check_numerical_error          = logical(1);
   inputs.read_bcs_from_initfile         = logical(0);
   inputs.use_capacitance_preconditioner = logical(1);
   inputs.use_deflation                  = logical(1);
   inputs.timesteps_between_writes       = 1;
   inputs.do_interfacial_averaging       = logical(0);

   % Write the input file.
   smpm_write_inputfile( configuration_file_name, inputs );

   % Write the tolerance file.
   fid = fopen( tolerance_file_name, 'w+' );
   fprintf( fid, '1e-10' );
   fclose( fid );

   % Write a dummy gold out file.
   fid = fopen( gold_file_name, 'w+' );
   fprintf( fid, 'This file is intentionally empty.' );
   fclose( fid );

return

function inputs_s = build_smpm_test_inputs( test_name )
% inputs_s = build_smpm_test_inputs( test_name )
%
% Creates an SMPM test case inputs structure that contains several default
% values and populates fields that are commonly derived from the test name.
% This is intended to be updated with configuration that should differ from
% the defaults (whatever they may be).
%
% Takes 1 value:
%
%   test_name - Name of the test.
%
% Returns 1 value:
%
%   inputs_s - Test case configuration structure that has been filled out.

   inputs_s.fname_runname = test_name;
   inputs_s.fname_init    = sprintf( '%s_init.h5',  test_name );
   inputs_s.fname_setup   = sprintf( '%s_setup', test_name );

   inputs_s.do_interfacial_averaging       = logical(0);
   inputs_s.check_null_error               = logical(0);
   inputs_s.check_numerical_error          = logical(0);
   inputs_s.read_bcs_from_initfile         = logical(0);
   inputs_s.use_capacitance_preconditioner = logical(0);
   inputs_s.use_deflation                  = logical(0);

   inputs_s.force_direct_solve             = logical(0);
   inputs_s.use_parallel_gmres             = logical(0);
   inputs_s.read_from_setupfile            = logical(0);
   inputs_s.write_to_setupfile             = logical(0);
   inputs_s.setup_and_stop                 = logical(0);

return
