// -*- adamgillard-cpp -*-
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
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

void FindNeighbours(std::vector<Boid>& flock)
{
    int num_boids = flock.size();
    double near_alert = 250.0; // Defines near distance
    double close_alert = 10.0; // Defines close distance
    double vision_ang = M_PI_2; // Defines angle over which a boid can see
    Vec3D disp_ij = Vec3D();

    // Perform Neighbour calculations over parallel i and sequential j (must perform over whole nxn matrix to prevent race conditions)
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int j = 0; j < num_boids; j++)
            {
                // Calculate distances between boids
                double dist = flock[i].pos.EuclidDist(flock[j].pos);

                if (dist != 0)    
                {
                    double mag = (flock[i].vel.Mag() * flock[j].vel.Mag());
                    double ang_ij = 0;
                    if (mag != 0)
                        ang_ij = acos(DotProd(flock[i].vel, flock[j].vel) / mag);

                    if (ang_ij < vision_ang)
                    {
                        // If close, append neighbour pos to boid close list 
                        if(dist < close_alert)
                        {   
                            // Calculate displacement 
                            disp_ij = flock[i].pos - flock[j].pos;
                            // Calculate mass-scaled position sum to use in calculating seperation force
                            flock[i].sum_pos_sep = flock[i].sum_pos_sep + (flock[j].mass * disp_ij);
                            // Add close-boid mass to close-mass sum and add 1 to number of detected close boids (for sanity check & debug)
                            flock[i].sum_cmass = flock[i].sum_cmass + flock[j].mass;
                            flock[i].num_close += 1;
                        }
                        // If not close but near, append neighbour pos to boid near list 
                        else if(dist < near_alert) 
                        {
                            // Calculate mass-scaled position sum to use for coherence force
                            flock[i].sum_pos_cohere = flock[i].sum_pos_cohere + ((1/flock[j].mass) * flock[j].pos);
                            // Calculate mass-scaled velocity sum to use for alignment force
                            flock[i].sum_vel_align = flock[i].sum_vel_align + (flock[j].mass * flock[j].vel);
                            // Add near-boid mass to close-mass sum and add 1 to number of detected near boids (for sanity check & debug)
                            flock[i].sum_nmass = flock[i].sum_nmass + flock[j].mass;
                            flock[i].num_near += 1;
                        }
                    }
                }
                
            } 
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
        CohereForce(flock[i]);
        SepForce(flock[i]);
        AlignForce(flock[i]);
        // Finally the bounding force from the walls
        WallForce(flock[i]);
        // Update the positions of the boids
        UpdatePos(flock[i]);
        // Reset the boid properties that vary each timestep
        Reset(flock[i]);

        // Make sure only one thread writes to file at a time using critical - may slow things down
        #pragma omp critical
        myfile << TIME << "," << iters << "," << flock[i].id << "," << flock[i].pos << "\n";
    } 
}

void ProgBar()
{
    // Start Progress Bar
    if (iters % 10 == 0)
    {
        std::cout << "[";
        for (int i = 0; i < 70; i++)
        {
            if (i < TIME/duration * 70) std::cout << "=";
            else if (i == TIME/duration * 70) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(TIME/duration * 100) << " %\r";
        std::cout.flush();
    }
    // End Progress Bar
}

int main(int argc, char* argv[])
{
    int argcheck = ArgChecks(argc, argv);
    if (argcheck)
        return 1;

    int threads = std::strtol(argv[1], nullptr, 0);
    int num_boids = std::strtol(argv[2], nullptr, 0);
    std::cout << "Program initialised for " << num_boids << " boids using " << threads << " OpenMP threads." << std::endl;

    omp_set_num_threads(threads); // Set OpenMP to use 10 threads
    auto tstart = std::chrono::high_resolution_clock::now();

    // Generate the flock in a vector on the heap
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<Boid> flock = GenFlock(num_boids, pos_ulim, pos_llim, vel_ulim, vel_llim);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
    std::cout << "Generated flock with " << flock.size() << " boids in " << time.count()/1e6 << " seconds.\n";

    std::ofstream myfile;
    myfile.open ("pos_dat.csv");
    myfile << "Time,Frame,ID,X,Y,Z,Type\n";

    Vec3D wind_force = RandWindForce();
    
    while (TIME < duration)
    {
        // Find all neighbours and store them in boid objects
        FindNeighbours(flock);
        // Complete force calculations, update positions and reset boids
        Simulate(flock, wind_force, myfile);
        // Evolve the wind force
        wind_force = WindEvo(wind_force);

        ProgBar();        

        TIME += dt;
        iters ++;
    }
    ProgBar();
    std::cout << std::endl;
    
    myfile.close();
    auto tend = std::chrono::high_resolution_clock::now();
    auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
    std::cout << "Completed Simulation for " << duration/dt << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";
    return 0;
}