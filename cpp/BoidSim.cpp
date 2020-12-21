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

void FindNeighbours(std::vector<Boid>& flock)
{
    int num_boids = flock.size();
    double near_alert = 300.0; // Defines near distance
    double close_alert = 10.0; // Defines close distance
    Vec3D sep_ij = Vec3D();

    for (int i = 0; i < num_boids-1; i++)
        for (int j = i+1; j < num_boids; j++)
            {
                // Calculate distances between boids
                double dist = flock[i].pos.EuclidDist(flock[j].pos);
                // If close, append neighbour pos to boid close list 
                if(dist < close_alert)
                {
                    sep_ij = flock[i].pos - flock[j].pos;
                    flock[i].sum_pos_sep = flock[i].sum_pos_sep + (flock[j].mass * sep_ij);
                    flock[j].sum_pos_sep = flock[j].sum_pos_sep + (flock[i].mass * (-1 * sep_ij));
                    flock[i].sum_cmass = flock[i].sum_cmass + flock[j].mass;
                    flock[j].sum_cmass = flock[j].sum_cmass + flock[i].mass;
                    flock[i].num_close += 1;
                    flock[j].num_close += 1;
                }
                // If not close but near, append neighbour pos to boid near list 
                else if(dist < near_alert) {
                    flock[i].sum_pos_cohere = flock[i].sum_pos_cohere + ((1/flock[j].mass) * flock[j].pos);
                    flock[j].sum_pos_cohere = flock[j].sum_pos_cohere + ((1/flock[i].mass) * flock[i].pos);
                    flock[i].sum_vel_align = flock[i].sum_vel_align + (flock[j].mass * flock[j].vel);
                    flock[j].sum_vel_align = flock[j].sum_vel_align + (flock[i].mass * flock[i].vel);
                    flock[i].sum_nmass = flock[i].sum_nmass + flock[j].mass;
                    flock[j].sum_nmass = flock[j].sum_nmass + flock[i].mass;
                    flock[i].num_near += 1;
                    flock[j].num_near += 1;
                }
            } 
}

void Simulate(std::vector<Boid>& flock, std::ofstream& myfile)
{
    for (int i = 0; i < num_boids; i++)
        {
            CohereForce(flock[i]);
            SepForce(flock[i]);
            AlignForce(flock[i]);
            WallForce(flock[i]);
            UpdatePos(flock[i]);
            Reset(flock[i]);
            myfile << TIME << "," << iters << "," << flock[i].id << "," << flock[i].pos << "\n";
        }
    ;
}


int main()
{
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
    
    
    while (TIME < duration)
    {
        // Find all neighbours and store them in boid objects
        
        FindNeighbours(flock);
        Simulate(flock, myfile);
        
        std::cout << "[";
        for (int i = 0; i < 70; i++)
        {
            if (i < TIME/duration * 70) std::cout << "=";
            else if (i == TIME/duration * 70) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(TIME/duration * 100) << " %\r";
        std::cout.flush();

        TIME += dt;
        iters ++;
    }
    std::cout << std::endl;
    
    myfile.close();
    auto tend = std::chrono::high_resolution_clock::now();
    auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
    std::cout << "Completed Simulation for " << duration/dt << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";
    return 0;
}