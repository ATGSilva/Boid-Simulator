/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Neighbours.cpp

Dependancy file to be used for building BoidSimOMP.cpp
Contains boid neighbour finding routines.

FUNCTION SIGNATURE - RETURN TYPE
    FindDists(vector<Boid>&) - vector<double>
    FindNeighbours(vector<Boid>&, vector<double>&) - void
    UpdateNeighboursFromBuffer(vector<Boid>&, vector<double>&) - void
*/

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


// Distance Finding ---------------------------------------------------
std::vector<double> FindDists(std::vector<Boid>& flock)
{
    /**
        Function to find NxN matrix of euclidean distances between all
        boids. Acts over upper triangular matrix due to symmetric
        nature of the problem.
        
        NOTES
        All instances of Timer are scope-based, hence the "random" use of
        scopes for single functions. See "Timing.h" for details.
    */
    int num_boids = flock.size();
    std::vector<double> distances(num_boids*num_boids, 0);

    #pragma omp parallel for
    for (int i = 0; i < num_boids-1; i++)
        for (int j = i+1; j < num_boids; j++)
        {
            double dist_ij = flock[i].pos.EuclidDist(flock[j].pos);
            distances[(i * num_boids) + j] = dist_ij;
            distances[(j * num_boids) + i] = dist_ij;
        }
            

    return distances;
}

void FindNeighbours(std::vector<Boid>& flock, std::vector<double>& dists)
{
    int num_boids = flock.size();

    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int j = 0; j < num_boids; j++)
        {
            double dist_ij = dists[(i *num_boids) + j];
                
            Vec3D rel_vec_ij = flock[j].pos - flock[i].pos; 
            double mag = (flock[i].vel.Mag() * rel_vec_ij.Mag());
            double ang_ij = 0;
            if (mag != 0)
                ang_ij = acos(DotProd(flock[i].vel, rel_vec_ij) / mag);

            if (ang_ij < VISION_FOV)
            {
                if (dist_ij < CLOSE_ALERT)
                    flock[i].close_list.push_back(flock[j].id);
                else if (dist_ij < NEAR_ALERT)
                    flock[i].near_list.push_back(flock[j].id);
            }

            if (dist_ij < BUFFER_ALERT)
                flock[i].buffer_list.push_back(flock[j].id);
        }
}

void UpdateNeighboursFromBuffer(std::vector<Boid>& flock, std::vector<double>& dists)
{
    int num_boids = flock.size();
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int id : flock[i].buffer_list)
        {
            double dist_iid = dists[(i *num_boids) + id];
            Vec3D rel_vec_iid = flock[id].pos - flock[i].pos; 
            double mag = (flock[i].vel.Mag() * rel_vec_iid.Mag());
            double ang_iid = 0;
            if (mag != 0)
                ang_iid = acos(DotProd(flock[i].vel, rel_vec_iid) / mag);

            if (ang_iid < VISION_FOV)
            {
                if (dist_iid < CLOSE_ALERT)
                    flock[i].close_list.push_back(flock[id].id);
                else if (dist_iid < NEAR_ALERT)
                    flock[i].near_list.push_back(flock[id].id);
            }
        }
}