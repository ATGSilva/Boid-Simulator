/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Simulate.cpp

Dependancy file to be used for building BoidSimMPIOMP.cpp
Contains boid neighbour finding routines and simulation wrapper.

FUNCTION SIGNATURE - RETURN TYPE
    FindDists(vector<Boid>&, int, int, int) - vector<double>
    FindNeighbours(vector<Boid>&, int, vector<double>&) - void
    Simulate(std::vector<Boid>&, int, int, std::vector<double>&, Vec3D, float)
*/

#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
#include "Settings.h"
#include "Forces.h"
#include "Boid.h"
#include "Vec3D.h"
#include "Simulate.h"

// Neighbour Finding --------------------------------------------------

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
    /**
        Function to find S*xN matrix of euclidean distances between all
        boids. Acts over whole matrix**.
        * S is the "split_by" matrix - each node only calculates over 
          an its equal share of the full boids, but must find the dists
          to all other boids in the flock.
        ** pseudo-matrix, the problem is dealt with as if it is acting
           over an SxN matrix but is in fact operating over a list of
           length S*N. A group of S elements is iterated over by i and 
           the number of N-element blocks into the list is given by j.
        
        OUTPUT
        List of length N*S of inter-boid distances.

        NOTES
        Nested loops are parallel over the inner loop only to prevent
        race condition.
    */

    int num_boids = flock.size();
    // Create vector of doubles and fill all S*N elements with zeros
    std::vector<double> distances(split_by*num_boids, 0);

    // Operate over whole matrix in parallel j (N elements)
    for (int i = 0; i < split_by; i++)
    {
        #pragma omp parallel for
        for (int j = 0; j < num_boids; j++)
        {
            // Find distance from ith to jth boid (always positive)
            // Map [i][j]th element to 1D list index
            distances[(i * num_boids) + j] = flock[flock_start+i].pos.EuclidDist(flock[j].pos);
        }
    }

    return distances;
}

float FovAngle(Boid& bi, Boid& bj)
{
    /**
        Function to find the angle of boid bi velocity vector to boid
        bj relative position vector. Range of arccos is 0 to pi hence
        this provides the absolute angle - this is fine as we don't
        care if it is a positive or negative angle as FOV is symmetric.
        
        INPUT
        bi - the boid we find the angle from.
        bj - the boid we find the angle to.

        OUTPUT
        float of angle value ranging between 0 and pi.
    */

    // Calculate the relative position vector from j to i
    Vec3D rel_vec_ij = bi.pos - bj.pos; 
    double mag = (bi.vel.Mag() * rel_vec_ij.Mag());
    float ang_ij = 0;
    if (mag != 0)
    {
        // Catch division by zero through if statement
        ang_ij = acos(DotProd(bi.vel, rel_vec_ij) / mag);
    }

    return ang_ij;
}

// Compiler instructions for progress bar control ---------------------
#define FOV VISION
#if FOV
#define FOVANGLE(Boid1, Boid2) FovAngle(Boid1, Boid2)
#else
#define FOVANGLE(Boid1, Boid2) 0
#endif
// If VISION == 1, make call to FOVANGLE() call FovAngle()
// Else, call to FOVANGLE() will always return 0


NeighbourLists FindNeighbours(std::vector<Boid>& flock, int i, std::vector<double>& dists)
{
    /**
        For each boid in the subflock, evaluate the distance to every
        other boid each check if it's within close or near ranges.
        If a boid is close (smallest radius category) then it is not
        included in the near list and vice versa - hence close boids
        only exert separation force, and near boids cohere and align.

        OUTPUTS
        Neighbour list instance containing near and close boids.
        
        NOTES
        VISION_FOV, CLOSE_ALERT, NEAR_ALERT defined in Settings.h.
        
        Must operate over all SxN boids since this is a non-symmetric
        problem due to the field-of-veiw calculation.
    */

    NeighbourLists lists;
    int num_boids = flock.size();

    // Outer loop is parallel and performed in Simulate() prior to
    // calling FindNeighbours()
    for (int j = 0; j < num_boids; j++)
    {   
        if (i != j)
        {
            // Read distance from ith to jth boid
            double dist_ij = dists[(i * num_boids) + j];
            // Calculate the angle from i velocity to j position
            float ang_ij = FOVANGLE(flock[i], flock[j]);

            // If j is within i field-of-view (i can see j) then append
            // j to neighbour list of i
            if (ang_ij < VISION_FOV)
            {
                if (dist_ij < CLOSE_ALERT)
                    lists.close_list.push_back(flock[j].id);
                else if (dist_ij < NEAR_ALERT)
                    lists.near_list.push_back(flock[j].id);
            }
        }
        
    }
    // If no neighbours at all (boid is alone), add to near i+1 boid 
    // so it moves towards the group
    if (lists.near_list.size() < 2)
        lists.near_list.push_back(flock[i+1].id);

    return lists;
}


// Simulation ---------------------------------------------------------

std::vector<Boid> Simulate(std::vector<Boid>& flock, int flock_start, int flock_end, std::vector<double>& distances, Vec3D wind_force)
{
    /**
        Runs all force calculations and position updates in parallel.
      
        @flock& memory reference to std::vector containing all boid objects.
        @flock_start first boid index to iterate through.
        @flock_end last boid index to iterate through.
        @distances& memory reference to inter-boid distance matrix.
        @wind_force& memory reference to wind_force variable.
    */

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

    return flock;
}