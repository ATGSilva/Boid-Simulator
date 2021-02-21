/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

BoidSimOMPMPI.cpp

Main file for running a flocking boid behavioural simulation.
This uses a hybrid MPI and OpenMP model for both multiprocessing and
multithreading.

FUNCTION SIGNATURE - RETURN TYPE
    Argchecks(int, char*) - int
    ProgBar() - void
    VecDataType() - MPI_Datatype
    BoidDataType() - MPI_Datatype
    Main(int, char*) - int

compile:
    mpicxx -fopenmp BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim -O3
run: 
    mpirun -np <NUM_PROCS> ./BoidSim <THREADS_PER_PROC> <NUM_BOIDS>
*/

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include "Simulate.h"
#include "TimingMPI.h"

// Global Variable Definitions -------------
float Time = 0;
int Iters = 0;
// Note: These are only global for this file.
//------------------------------------------

// Utilities ----------------------------------------------------------

int ArgChecks(int argc, char* argv[])
{
    /**
        Performs checks on all execution input variables to ensure the program
        will not crash due to misuse.
      
        @argc The number of input parameters supplied on execution.
        @argv pointer to start element of array containing supplied input 
        values.
    */

    if (argc != 3)
    {
        std::cerr << "Usage: ./BoidSimOMP <Integer NUM_THREADS> <Integer NUM_BOIDS>" << std::endl;
        return 1;
    }

    int threads = std::strtol(argv[1], nullptr, 0);
    if (threads < 1 || threads > omp_get_max_threads())
    {
        std::cerr << "<NUM_THREADS> must be between 1 and " << omp_get_max_threads() << std::endl;
        return 1;
    }
    int num_boids = std::strtol(argv[2], nullptr, 0);
    if (num_boids < 2)
    {
        std::cerr << "<NUM_BOIDS> must be a positive integer more than 1." << std::endl;
        return 1;
    }
    else return 0;
}

// Compiler instructions for progress bar control
#define PROGRESSUPDATE PROGRESS
#if PROGRESSUPDATE
#define PROGBAR() ProgBar()
#else
#define PROGBAR()
#endif

void ProgBar()
{
    // Start Progress Bar
    std::cout << "[";
    for (int i = 0; i < 70; i++)
    {
        if (i < Time/DURATION * 70) std::cout << "=";
        else if (i == Time/DURATION * 70) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(Time/DURATION * 100) << " %\r";
    std::cout.flush();
    // End Progress Bar
}

// MPI Datatypes ------------------------------------------------------

MPI_Datatype VecDataType()
{
    Vec3D vec = Vec3D(1., 2., 3.);
    MPI_Datatype VECTYPE;
    int vlens[3] = {1,1,1};
    MPI_Aint vdisps[3] = {0};
    MPI_Aint vst_addr, vaddr;
    MPI_Get_address(&vec.x, &vst_addr);
    MPI_Get_address(&vec.y, &vaddr);
    vdisps[1] = vaddr - vst_addr;
    MPI_Get_address(&vec.z, &vaddr);
    vdisps[2] = vaddr - vst_addr;
    MPI_Datatype vtypes[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(3, vlens, vdisps, vtypes, &VECTYPE);
    MPI_Type_commit(&VECTYPE);

    return VECTYPE;
}

MPI_Datatype BoidDataType(MPI_Datatype VECTYPE)
{
    Boid boid = Boid(Vec3D(1.0f,2.0f,3.0f), Vec3D(4.0f,5.0f,6.0f), 0);
    MPI_Datatype BOIDTYPE;
    int blens[4] = {1, 1, 1, 1};
    MPI_Aint bdisps[4] = {0};
    MPI_Aint bst_addr, baddr;
    MPI_Get_address(&boid.pos, &bst_addr);
    MPI_Get_address(&boid.vel, &baddr);
    bdisps[1] = baddr - bst_addr;
    MPI_Get_address(&boid.mass, &baddr);
    bdisps[2] = baddr - bst_addr;
    MPI_Get_address(&boid.id, &baddr);
    bdisps[3] = baddr - bst_addr;
    MPI_Datatype btypes[4] = {VECTYPE, VECTYPE, MPI_FLOAT, MPI_INT};
    MPI_Type_create_struct(4, blens, bdisps, btypes, &BOIDTYPE);
    MPI_Type_commit(&BOIDTYPE);

    return BOIDTYPE;
}

// Compiler instructions for readout file control ---------------------
#define PROFILING BOID_READOUT
#if PROFILING
#define START_PSESSION(results, num_boids, numnodes, threads) BeginResultSession(results, num_boids, numnodes, threads)
#define PROFILE_READOUT(results, boid) WriteResults(results, Time, Iters, boid)
#define END_PSESSION(results) EndWriteSession(results)
#else
#define START_PSESSION(results, num_boids, numnodes, threads)
#define PROFILE_READOUT(results, boid)
#define END_PSESSION(results)
#endif

#define BENCHMARK TIMING
#if BENCHMARK
#define START_TSESSION(timingfile, num_boids, numnodes, threads) BeginTimingSession(timingfile, num_boids, numnodes, threads)
#define TIMER(name, timingfile) Timer timer(name, timingfile, Iters)
#define END_TSESSION(timingfile) EndWriteSession(timingfile)
#else
#define START_TSESSION(timingfile, num_boids, numnodes, threads)
#define TIMER(name, timingfile)
#define END_TSESSION(timingfile)
#endif


// Main Functionality -------------------------------------------------

int main(int argc, char* argv[])
{
    /**
        Program to simulate the flocking behaviour of boids over N-bodies.
        Parallelisation by MPI and OpenMP.
        
        NOTES
        All instances of Timer are scope-based, hence the "random" use of
        scopes for single functions. See "Timing.h" for details.
    */

   // INPUT CHECKS ----------------------------------------------------
    int argcheck = ArgChecks(argc, argv);
    if (argcheck)
        return 1;

    // Variable Definitions & Declarations ----------------------------
    int threads = std::strtol(argv[1], nullptr, 0);
    int num_boids = std::strtol(argv[2], nullptr, 0);
    int num_timesteps = DURATION / DT;
    Vec3D wind_force;
    std::vector<double> dist_matrix;
    std::vector<Boid> flock;
    int flock_start, flock_end, split_by, remaining, numtasks, numworkers, taskid;

    // OpenMP Setup ---------------------------------------------------
    omp_set_num_threads(threads); // Set OpenMP to use input no of threads

    // MPI Setup ------------------------------------------------------
    
    MPI_Init(&argc, &argv);                     // Initialise the MPI environment
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);   // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);     // Get the rank of the process
    numworkers = numtasks - 1;                  // Calculate number of workers
    MPI_Status status;                          // Define MPI status variable

    std::vector<int> BEGIN_TAGS(num_timesteps);                         // 
    std::vector<int> WIND_TAGS(num_timesteps);                          // Set up tag lists
    std::vector<int> DONE_TAGS(num_timesteps);                          //
    std::iota(BEGIN_TAGS.begin(), BEGIN_TAGS.end(), 0);                 // 
    std::iota(WIND_TAGS.begin(), WIND_TAGS.end(), num_timesteps+1);     // Lists of incrementing integers 
    std::iota(DONE_TAGS.begin(), DONE_TAGS.end(), 2*(num_timesteps+1)); //

    MPI_Datatype VECTYPE = VecDataType();           // Set up MPI Datatypes
    MPI_Datatype BOIDTYPE = BoidDataType(VECTYPE);  //

    
    split_by = num_boids / numtasks;    // Set up factors to determine how many 
    remaining = num_boids % numtasks;   // boids are dealt with on each core
    if (remaining != 0)
    {
        // Force the number of boids to be a multiple of the number of cores used for easy splitting 
        num_boids += (numtasks - remaining);
        if (taskid == DIRECTOR)
        {
            std::cout << "Adding " << (numtasks - remaining) << " boids to allow even flock splitting." << "\n";
        }
        split_by = num_boids / numtasks; // Recalculate split_by
    }

    // Other Setup ----------------------------------------------------
    std::ofstream results;          // Set up output files
    std::ofstream timingfile;       //
    auto tstart = std::chrono::high_resolution_clock::now();

    // Director Specific Setup ----------------------------------------
    if (taskid == DIRECTOR)
    {
        if (numworkers > MAXWORKER || numworkers < MINWORKER)
        {
            std::cout << "The number of MPI Processors must be between " << MINWORKER+1 << " and " << MAXWORKER+1 << "\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::cout << "Program initialised for " << num_boids << " boids using " << (threads*numtasks) << " OpenMP threads across " << numtasks << " MPI nodes.\n";

        // Generate the flock in a vector on the heap
        flock.reserve(num_boids);
        flock = GenFlock(num_boids, POS_ULIM, POS_LLIM, VEL_ULIM, VEL_LLIM);

        // Open data files and give headers
        START_PSESSION(results, num_boids, numtasks, threads);
        START_TSESSION(timingfile, num_boids, numtasks, threads);

        // Generate a random wind force
        wind_force = RandWindForce();
    }

    // Syncronise and run through all timesteps -----------------------
    MPI_Barrier(MPI_COMM_WORLD);
    while (Time < DURATION)
    {
        TIMER("Timestep", timingfile);

        if (taskid == DIRECTOR)
        {
            // BEGIN DIRECTOR CODE ------------------------------------

            // Send full flock and wind force to all workers
            for (int w_id = 1; w_id < numworkers+1; w_id++)
            {
                MPI_Send(&flock.front(), num_boids, BOIDTYPE, w_id, BEGIN_TAGS[Iters], MPI_COMM_WORLD);
                MPI_Send(&wind_force, 1, VECTYPE, w_id, WIND_TAGS[Iters], MPI_COMM_WORLD);
            }

            // Simulate Boid evo for first <split_by> Boids
            flock_start = 0;
            flock_end = split_by;
            dist_matrix.resize(split_by * num_boids);
            {
                TIMER("Distances", timingfile);
                dist_matrix = FindDists(flock, flock_start, flock_end, split_by);
            }
            {
                TIMER("Simulate", timingfile);
                flock = Simulate(flock, flock_start, flock_end, dist_matrix, wind_force, Time);
            }
            
            // Recieve simulated subflocks back from workers into flock
            for (int w_id = 1; w_id < numworkers+1; w_id++)
            {
                int offset_start = (w_id) * split_by;
                MPI_Recv(&flock[offset_start], split_by, BOIDTYPE, w_id, DONE_TAGS[Iters], MPI_COMM_WORLD, &status);
            }

            // Write data to file
            for (int num = 0; num < num_boids; num++)
                PROFILE_READOUT(results, flock[num]);

            // Evolve the wind force
            wind_force = WindEvo(wind_force);
			
            // Update Progress bar every 20 iterations
			if (Iters % 20 == 0)
			{
				PROGBAR();   
			}

            // END DIRECTOR CODE --------------------------------------
        }
        else
        {
            // BEGIN WORKER CODE --------------------------------------
            
            // Set up variables ready for recieve
            flock.resize(num_boids);
            Vec3D wind_force;
            flock_start = (taskid) * split_by;
            flock_end = flock_start + split_by;

            // Recieve flock and wind force from the Director
            MPI_Recv(&flock.front(), num_boids, BOIDTYPE, DIRECTOR, BEGIN_TAGS[Iters], MPI_COMM_WORLD, &status);
            MPI_Recv(&wind_force, 1, VECTYPE, DIRECTOR, WIND_TAGS[Iters], MPI_COMM_WORLD, &status);

            // Perform calculation over a subset of the flock
            dist_matrix.resize(split_by * num_boids);
            dist_matrix = FindDists(flock, flock_start, flock_end, split_by);
            flock = Simulate(flock, flock_start, flock_end, dist_matrix, wind_force, Time);

            // Send the simulated subflock back to Director
            MPI_Send(&flock[flock_start], split_by, BOIDTYPE, DIRECTOR, DONE_TAGS[Iters], MPI_COMM_WORLD);

            // Clear flock to ensure correct recive on next timestep
            flock.clear();

            // END WORKER CODE ----------------------------------------
        }

        // Syncronise at end of timestep
        MPI_Barrier(MPI_COMM_WORLD);
        Time += DT; // Update time by timestep
        Iters ++; // Add 1 to the iteration value
    }
    
    
    // Clean up ------------------------------------------------------- 

    // Release the memory for custom datatypes
    MPI_Type_free(&BOIDTYPE);
    MPI_Type_free(&VECTYPE);
    // Finalize the MPI environment.
    MPI_Finalize();

    if (taskid == DIRECTOR)
    {
        PROGBAR();
        END_PSESSION(results);
        END_TSESSION(timingfile);
        auto tend = std::chrono::high_resolution_clock::now();
        auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
		WriteSettings(num_boids, numtasks, threads);
        std::cout << "\nCompleted Simulation for " << DURATION/DT << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";
    }

    return 0;
}