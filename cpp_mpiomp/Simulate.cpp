#include "Settings.h"
#include "Forces.h"
#include "Boid.h"
#include "Vec3D.h"
#include "Simulate.h"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>

void ProgBar(float& time)
{
    // Start Progress Bar
    std::cout << "[";
    for (int i = 0; i < 70; i++)
    {
        if (i < time/duration * 70) std::cout << "=";
        else if (i == time/duration * 70) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(time/duration * 100) << " %\r";
    std::cout.flush();
    // End Progress Bar
}

struct NeighbourLists
{
    NeighbourLists() {
        std::vector<int> close_list(1);
        std::vector<int> near_list(1);
    };

    std::vector<int> close_list;
    std::vector<int> near_list;
};

std::vector<double> FindDists(std::vector<Boid>& flock, int flock_start, int flock_end, int split_by)
{
    int num_boids = flock.size();
    std::vector<double> distances;
    distances.reserve(num_boids * split_by);

    for (int i = 0; i < split_by; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < num_boids; j++)
            distances[(i *num_boids) + j] = flock[flock_start+i].pos.EuclidDist(flock[j].pos);
    }

    return distances;
}

NeighbourLists FindNeighbours(std::vector<Boid>& flock, int i, std::vector<double>& dists)
{
    NeighbourLists lists;
    int num_boids = flock.size();

    for (int j = 0; j < num_boids; j++)
    {
        double dist_ij = dists[(i * num_boids) + j];

        Vec3D rel_vec_ij = flock[j].pos - flock[i].pos; 
        double mag = (flock[i].vel.Mag() * rel_vec_ij.Mag());
        double ang_ij = 0;
        if (mag != 0)
            ang_ij = acos(DotProd(flock[i].vel, rel_vec_ij) / mag);

        if (ang_ij < vision_ang)
        {
            if (dist_ij < close_alert)
                lists.close_list.push_back(flock[j].id);
            else if (dist_ij < near_alert)
                lists.near_list.push_back(flock[j].id);
        }
    }
    // If no neighbours at all (boid is alone), add to near i+1 boid so it moves towards the group
    if (lists.near_list.size() < 1)
        lists.near_list.push_back(flock[i+1].id);

    return lists;
}

std::vector<Boid> Simulate(std::vector<Boid>& flock, int flock_start, int flock_end, std::vector<double>& distances, Vec3D wind_force, float time)
{
    #pragma omp parallel for
    for (int i = flock_start; i < flock_end; i++)
    {
        // Determine neighbours
        NeighbourLists neighbour_lists = FindNeighbours(flock, i-flock_start, distances);
        // Calculate forces
        Forces force_list;
        force_list.cohere_force = CohereForce(flock, flock[i], neighbour_lists.near_list);
        force_list.sep_force = SepForce(flock, flock[i], neighbour_lists.close_list);
        force_list.align_force = AlignForce(flock, flock[i], neighbour_lists.near_list);
        force_list.wind_force = wind_force; 
        force_list.wall_force = WallForce(flock[i]);
        // Update boid position
        UpdatePos(flock[i], force_list);
    }
    ProgBar(time);

    return flock;
}