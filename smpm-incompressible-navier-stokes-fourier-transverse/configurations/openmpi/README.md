# Overview

This directory contains OpenMPI-specific configuration files.  These aren't
necessary for execution on all systems, though if problems arise when using the
testing frameworks (benchmarks or validation), these should be used as a
workaround until a newer version of OpenMPI (no older than 2.0) can be
installed.

# Installation

The OpenMPI MCA parameter file can be used as is via a run-time argument like
so:

``` shell
$ mpirun -np 16 -am path/to/mca-params.conf ./smpm_incompressible_navier_stokes path/to/run_in path/to/run_out path/to/run_restart
```

It may also be installed be in a more global manner, either system-wide (which
requires `root` and is expected to be infrequently used) or user-specific:

``` shell
$ mkdir ~/.openmpi
$ cp mca-params.conf ~/.openmpi
```
