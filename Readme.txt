=========================================================================================
  ____  _____    ____        _     _    _____ _                 _       _             
 |___ \|  __ \  |  _ \      (_)   | |  / ____(_)               | |     | |            
   __) | |  | | | |_) | ___  _  __| | | (___  _ _ __ ___  _   _| | __ _| |_ ___  _ __ 
  |__ <| |  | | |  _ < / _ \| |/ _` |  \___ \| | '_ ` _ \| | | | |/ _` | __/ _ \| '__|
  ___) | |__| | | |_) | (_) | | (_| |  ____) | | | | | | | |_| | | (_| | || (_) | |   
 |____/|_____/  |____/ \___/|_|\__,_| |_____/|_|_| |_| |_|\__,_|_|\__,_|\__\___/|_|   

=========================================================================================

Author: Adam Gillard
Unit: Advanced Computational Physics - University of Bristol
Contact: ag17009@bristol.ac.uk
Alt Contact: adamgillard1@hotmail.co.uk


CONTENTS
--------
 * Description
 * Versions
 * Requirements
 * Configuration
 * Execution
 * Animation & Dashboard Plotting (IMPORTANT READ)
 * Clean


DESCRIPTION
-----------

A program to simulate swarming and flocking of bird-oids (boids) in 3D space, based on 
the rules created by Craig Reynolds (https://en.wikipedia.org/wiki/Boids).


VERSIONS
--------

 1) Multithreading with OpenMP
 2) Hybrid multiprocessing & multithreading with MPI & OpenMP


REQUIREMENTS
------------

These programs have been tested on: 
 - "Ubuntu 20.04.2 LTS" though WSL2 (AMD Ryzen 5 2600)
 - Blue Crystal Phase 3

This project uses a makefile system for building, running and plotting.

For OpenMP version, makefile build is set up to use gcc, and has been tested with:
 - GCC v9.3.0
 - GCC v4.8.5

For MPI version, makefile build is set up to use mpicxx, and has been tested with:
 - MPICXX -> GCC v7.3.0
 - MPICXX -> GCC v4.8.5
 
Blue Crystal Modules:
 - languages/gcc-4.8.5
 - openmpi/gcc/64/1.6.5

WARNING: Ensure these are the ONLY LOADED MODULES to prevent conflicts.

WARNING: Plotting and Animations will not work on Blue Crystal! Python modules
		 cause disastrous conflicts with mpi program execution! 


CONFIGURATION
-------------

Before running a program, you may wish to edit the configuration file to change the
behaviour of the boids. This can be done by editing the file:
 - "PATH/src/<omp/mpi>/Settings.h"
where <omp/mpi> is a choice depending on which version of the program you wish to run.

Current settings can be checked on linux command line using
 `cat src/<omp/mpi>/Settings.h`
and is loaded with default values that produce almost-perfect flocking.

Settings can be changed in command line similarly to above using your editor or choice.


EXECUTION 
---------

When executing on a personal computer, you are free to use all features of the makefile.
When executing on Blue Crystal, ONLY use the build functionality of makefile!
 - Runscripts are provided that can be edited to suit run configurations.

To build programs:
 - Build OpenMP program
   - `make build_omp`
   or
   - `cd src/omp`
   and
   - `g++ <FLAGS> BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp Vec3D.cpp 
      -o ../../BoidSimOMP`
 - Build MPI Program:
   - `make build_mpi`
   or
   - `cd src/mpi`
   and
   - `mpicxx <FLAGS> BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp 
      -o ../../BoidSimMPI`

Note - manual build options can be used to add machine-specific flags.

To run programs (from main directory):
`make run_omp (NUM_BOIDS=<desired> OMP_THR=<number of threads>)`
`make run_mpi (NUM_BOIDS=<desired> NP=<number of nodes> PPN=<threads per node>)`
or
`./BoidSimOMP <number of threads> <number of boids>` == OpenMP
`mpirun -np <number of nodes> ./BoidSimMPI <threads per node> <number of boids>` == MPI


ANIMATION & DASHBOARD PLOTTING
------------------------------

Animations are 3D visualisations of the boids over the entire run duration.
Dashboards include the above alongside animated average momentum and flocking
coefficient plots.

WARNING: These plot every single timestep as a frame. If you are using short timesteps
(<0.5) you will need to change the `skip_factor` and/or `fps` at the top of the
"utilities/<omp/mpi>/plot_*.py" files for reasonable generation time and playback speed.

Animations and dashboards read in the "boid_data/settings.cfg" file produced by the most
recent execution! It uses this information to read in the correct boid data files, label
axes and produce flocking coefficients. If you wish to keep settings files for later
plotting, create a copy of the file for later use - when you wish to use a file, make
sure it is "boid_data/settings.cfg".

If you wish to plot data produced on Blue Crystal, you must copy the following files to
your local PC alongside the utilities in the same directory structure:
 - For OpenMP files:
   - boid_data/B<num_boids>Thr<num_threads>.csv
   - boid_data/settings.cfg (corresponding to the correct execution run)
 - For MPI files:
   - boid_data/B<num_boids>Nodes<num_nodes>Thr<threads_per_node>.csv
   - boid_data/settings.cfg (corresponding to the correct execution run)
These must be copied into your local "boid_data" directory, which is created when you
build a program via the makefile.

To create animations of the boid simulation:
`make animation_<omp/mpi>`

To create a dashboard animation (3D plot, average momentum, flocking coefficients):
`make dashboard_<omp/mpi>`
   

CLEAN
-----

The makefile option `make clean` will remove:
 - The binary executables
 - ALL contents of "timing/" including the directory
 - ALL contents of "boid_data/" including the directory
 - ALL contents of "animations/" including the directory


