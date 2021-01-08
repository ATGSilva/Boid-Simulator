// -*- adamgillard-cpp -*-
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
//#include <mpi.h>
#include "Neighbours.h"
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"

float TIME = 0;
int iters = 0;

int ArgChecks(int argc, char* argv[])
{
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

void Simulate(std::vector<Boid>& flock, Vec3D& wind_force, std::ofstream& myfile)
{   
    int num_boids = flock.size();
    // Calculate forces and movement updates on every boid in the flock parallel-wise
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
    {
        // Wind flock affects the entire flock by the same amount
        flock[i].wind_force = wind_force;
        // Other forces are dependent on the neighbours
        CohereForce(flock, flock[i]);
        SepForce(flock, flock[i]);
        AlignForce(flock, flock[i]);
        // Finally the bounding force from the walls
        WallForce(flock[i]);
        // Update the positions of the boids
        UpdatePos(flock[i]);
    } 
    // To prevent race conditions and limit use of critical (slow), reloop through all sequentially to write to file and ensure proper reset
    for (int i = 0; i < num_boids; i++)
    {
        // Reset the boid properties that vary each timestep
        Reset(flock[i], iters);
        myfile << TIME << "," << iters << "," << flock[i].id << "," << flock[i].pos << "\n";
    }
}

void ProgBar()
{
    // Start Progress Bar
    std::cout << "[";
    for (int i = 0; i < 70; i++)
    {
        if (i < TIME/duration * 70) std::cout << "=";
        else if (i == TIME/duration * 70) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(TIME/duration * 100) << " %\r";
    std::cout.flush();
    // End Progress Bar
}


int main(int argc, char* argv[])
{
    int argcheck = ArgChecks(argc, argv);
    if (argcheck)
        return 1;

    //MPI_Init(NULL, NULL);

    // Get the number of processes
    //int world_size;
    //MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    //int world_rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    //char processor_name[MPI_MAX_PROCESSOR_NAME];
    //int name_len;
    //MPI_Get_processor_name(processor_name, &name_len);

    // Finalize the MPI environment.
    //MPI_Finalize();


    int threads = std::strtol(argv[1], nullptr, 0);
    int num_boids = std::strtol(argv[2], nullptr, 0);
    std::cout << "Program initialised for " << num_boids << " boids using " << threads << " OpenMP threads." << std::endl;

    omp_set_num_threads(threads); // Set OpenMP to use input no of threads
    auto tstart = std::chrono::high_resolution_clock::now();

    // Generate the flock in a vector on the heap
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<Boid> flock = GenFlock(num_boids, pos_ulim, pos_llim, vel_ulim, vel_llim);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
    std::cout << "Generated flock with " << flock.size() << " boids in " << time.count()/1e6 << " seconds.\n";

    std::ofstream myfile;
    myfile.open ("pos_dat.csv");
    myfile << "Time,Frame,ID,X,Y,Z\n";

    std::ofstream timingfile;
    timingfile.open ("timing.csv");
    timingfile << "Frame,Neighbour,Simulation\n";

    Vec3D wind_force = RandWindForce();
    std::vector<double> distances;

    while (TIME < duration)
    {
        // Find all neighbours and store them in boid objects
        auto t1 = std::chrono::high_resolution_clock::now();
        if (iters % buffer_for == 0)
        {
            
            distances = FindDists(flock);
            FindNeighbours(flock, distances);
            
        }
        else
        {
            UpdateNeighboursFromBuffer(flock, distances);
        }
        auto t2 = std::chrono::high_resolution_clock::now();    

        // Complete force calculations, update positions and reset boids
        Simulate(flock, wind_force, myfile);
        auto t3 = std::chrono::high_resolution_clock::now(); 
        // Evolve the wind force
        wind_force = WindEvo(wind_force);

        auto time_neighb = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
        auto time_sim = std::chrono::duration_cast<std::chrono::microseconds>(t3-t2);
        timingfile << iters << "," << time_neighb.count()/1e6 << "," << time_sim.count()/1e6 << "\n";
        
        if (iters % 10 == 0)
            ProgBar();        

        TIME += dt;
        iters ++;
    }
    ProgBar();
    std::cout << std::endl;
    
    myfile.close();
    timingfile.close();
    auto tend = std::chrono::high_resolution_clock::now();
    auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
    std::cout << "Completed Simulation for " << duration/dt << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";
    return 0;
}