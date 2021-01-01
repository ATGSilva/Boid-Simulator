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

std::vector<double> FindDists(std::vector<Boid>& flock)
{
    int num_boids = flock.size();
    std::vector<double> distances;
    distances.reserve(num_boids * num_boids);

    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int j = 0; j < num_boids; j++)
            distances[(i *num_boids) + j] = flock[i].pos.EuclidDist(flock[j].pos);

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

            if (ang_ij < vision_ang)
            {
                if (dist_ij < close_alert)
                    flock[i].close_list.push_back(flock[j].id);
                else if (dist_ij < near_alert)
                    flock[i].near_list.push_back(flock[j].id);
            }

            if (dist_ij < buffer_alert)
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

            if (ang_iid < vision_ang)
            {
                if (dist_iid < close_alert)
                    flock[i].close_list.push_back(flock[id].id);
                else if (dist_iid < near_alert)
                    flock[i].near_list.push_back(flock[id].id);
            }
        }
}