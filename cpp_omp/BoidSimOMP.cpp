/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

BoidSimOMP.cpp

Main file for running a flocking boid behavioural simulation.

FUNCTION SIGNATURE - RETURN TYPE
    Argchecks(int, char*) - int
    ProgBar() - void
    Simulate(vector<Boid>&, Vec3D&, ofstream&) - void
    Main(int, char*) - int
*/

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <string>
#include "Neighbours.h"
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include "Timing.h"

// Global Variable Definitions -------------
float Time = 0;
int Iters = 0;
// Note: These are only global for this file.
//------------------------------------------


// Utilities ------------------------------------------------------------------------------------------------------

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


// Compiler instructions for readout file control ---------------------
#define PROFILING 1
#if PROFILING
#define START_PSESSION(results) BeginResultSession(results)
#define PROFILE_READOUT(results, boid) WriteResults(results, Time, Iters, boid)
#define END_PSESSION(results) EndWriteSession(results)
#else
#define START_PSESSION(results)
#define PROFILE_READOUT(results, boid)
#define END_PSESSION(results)
#endif

#define BENCHMARK 1
#if BENCHMARK
#define START_TSESSION(timingfile, num_boids, threads) BeginTimingSession(timingfile, num_boids, threads)
#define TIMER(name, timingfile) Timer timer(name, timingfile, Iters)
#define END_TSESSION(timingfile) EndWriteSession(timingfile)
#else
#define START_TSESSION(timingfile, num_boids, threads)
#define TIMER(name, timingfile)
#define END_TSESSION(timingfile)
#endif


// Main Functionality -------------------------------------------------

void Simulate(std::vector<Boid>& flock, Vec3D& wind_force, std::ofstream& results) 
{
    /**
        Runs all force calculations and position updates in parallel. Resets boid
        properties that vary on each timestep.
      
        @flock& memory reference to std::vector containing all boid objects.
        @wind_force& memory reference to wind_force variable.
        @results& memory reference to file which records time, frame and boid 
        position.
    */
    int num_boids = flock.size();
    
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
    {
        // Loop over every boid i in the flock, and calculate forces and update position
        flock[i].wind_force = wind_force; // Wind flock affects the entire flock by the same amount
        
        CohereForce(flock, flock[i]);   //
        SepForce(flock, flock[i]);      // Calculate core behavioural forces
        AlignForce(flock, flock[i]);    //
        WallForce(flock[i]);            // Finally the bounding force from the walls
        UpdatePos(flock[i]);            // Update the positions of the boids
    }

    // To prevent race conditions and limit use of critical (slow)
    // reloop through all sequentially to write to file and ensure proper reset
    for (int i = 0; i < num_boids; i++)
    {
        Reset(flock[i], Iters); // Reset the boid properties that vary each timestep
        // Write "current time, iteration number, boid ID number, boid position, boid velocity" to results file
        PROFILE_READOUT(results, flock[i]);
        //results << Time << "," << Iters << "," << flock[i].id << "," << flock[i].mass << "," << flock[i].pos << "," << flock[i].vel <<"\n";
    }
}

int main(int argc, char* argv[]) 
{
    /**
        Program to simulate the flocking behaviour of boids over N-bodies.
        
        NOTES
        All instances of Timer are scope-based, hence the "random" use of
        scopes for single functions. See "Timing.h" for details.
    */

    // INPUT CHECKS -----------------------------------------
    int argcheck = ArgChecks(argc, argv);
    if (argcheck) 
        return 1;

    // Variable Definitions & Declarations ------------------
    const int threads = std::strtol(argv[1], nullptr, 0);
    const int num_boids = std::strtol(argv[2], nullptr, 0);
    std::vector<Boid> flock;
    Vec3D wind_force;
    std::vector<double> distances;
    std::ofstream results, timingfile;

    // Other Setup ------------------------------------------
    START_PSESSION(results);
    START_TSESSION(timingfile, num_boids, threads);
    //BeginTimingSession(timingfile, num_boids, threads);                     // Set up timing file ready to be written to
    //BeginResultSession(results);                                            // Set up results file ready to be written to
    omp_set_num_threads(threads);                                           // Set OpenMP to use input no of threads
    flock = GenFlock(num_boids, POS_ULIM, POS_LLIM, VEL_ULIM, VEL_LLIM);    // Generate a flock of boids
    wind_force = RandWindForce();                                           // Generate a random initial wind force

    auto tstart = std::chrono::high_resolution_clock::now();
    std::cout << "Program initialised for " << num_boids << " boids using " << threads << " OpenMP threads.\n";

    // Run through all timesteps ---------------------------------------
    while (Time < DURATION)
    {
        Timer timer("Timestep", timingfile, Iters);

        // Find all neighbours and store them in boid objects
        if (Iters % BUFFER_FOR == 0) // Don't use buffer
        {    
            { // Timing Scope
            TIMER("Distances", timingfile);
            //Timer timer("Distances", timingfile, Iters);
            distances = FindDists(flock);
            }
            TIMER("Neighb", timingfile);
            //Timer timer("Neighb", timingfile, Iters);
            FindNeighbours(flock, distances);
        }
        else // Use buffer
        {
            TIMER("NeighbBuffer", timingfile);
            //Timer timer("NeighbBuffer", timingfile, Iters);
            UpdateNeighboursFromBuffer(flock, distances);
        }

        // Complete force calculations, update positions and reset boids
        { // Timing Scope
        TIMER("Simulate", timingfile);
        //Timer timer("Simulate", timingfile, Iters);
        Simulate(flock, wind_force, results);
        }
        
        // Update Progress bar every 20 iterations
        if (Iters % 20 == 0)
        {
            ProgBar();   
        }

        wind_force = WindEvo(wind_force); // Evolve the wind force
        Time += DT; // Update time by timestep
        Iters ++; // Add 1 to the iteration value
    }
    
    // Perform IO shutdown ----------------------------------
    ProgBar();
    END_PSESSION(results);
    END_TSESSION(timingfile);
    //EndWriteSession(timingfile);
    //EndWriteSession(results);

    auto tend = std::chrono::high_resolution_clock::now();
    auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
    std::cout << "\nCompleted Simulation for " << DURATION/DT << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";

    return 0;
}