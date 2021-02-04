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
        boids. Acts over upper triangular matrix* due to symmetric
        nature of the problem reducing complexity to O(nlogn).
        * pseudo-matrix, the problem is dealt with as if it is acting
          over an NxN matrix but is in fact operating over a list of
          length N^2. A group of N elements is iterated over by i and 
          the number of N-element blocks into the list is given by j.
        
        OUTPUT
        List of length N^2 of inter-boid distances.

        NOTES
        Nested loops are parallel over the outer loop only to prevent
        race condition.
    */

    int num_boids = flock.size();
    // Create vector of doubles and fill all N*N elements with zeros
    std::vector<double> distances(num_boids*num_boids, 0);

    // Operate over upper triagle - not including diagonal
    #pragma omp parallel for
    for (int i = 0; i < num_boids-1; i++)
        for (int j = i+1; j < num_boids; j++)
        {
            // Find distance from ith to jth boid (always positive)
            double dist_ij = flock[i].pos.EuclidDist(flock[j].pos);
            // Map [i][j]th and [j][i]th elements to 1D list index
            distances[(i * num_boids) + j] = dist_ij;
            distances[(j * num_boids) + i] = dist_ij;
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

void FindNeighbours(std::vector<Boid>& flock, std::vector<double>& dists)
{
    /**
        For each boid in the flock, evaluate the distance to every
        other boid each check if it's within close, near, or buffer.
        If a boid is close (smallest radius category) then it is not
        included in the near list and vice versa - hence close boids
        only exert separation force, and near boids cohere and align.
        Buffer lists include all buffer, near and close boids.

        OUTPUTS
        Neighbour lists are stored within the Boid object as a class
        property.
        
        NOTES
        VISION_FOV, CLOSE_ALERT, NEAR_ALERT, BUFFER_ALERT defined in
        Settings.h.
        Must operate over all NxN boids now since this is a non-
        symmetric problem due to the field-of-veiw calculation.
    */

    int num_boids = flock.size();

    // Operate over all NxN boids - non-symmetric problem
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int j = 0; j < num_boids; j++)
        {
            // Don't consider i = j case (boids i and j are the same)
            if (i != j)
            {
                // Read distance from ith to jth boid
                double dist_ij = dists[(i *num_boids) + j];

                // Calculate the angle from i velocity to j position
                float ang_ij = FovAngle(flock[i], flock[j]);

                // If j is within i field-of-view (i can see j) then append
                // j to corresponding neighbour list of i
                if (ang_ij < VISION_FOV)
                {
                    if (dist_ij < CLOSE_ALERT)
                        flock[i].close_list.push_back(flock[j].id);
                    else if (dist_ij < NEAR_ALERT)
                        flock[i].near_list.push_back(flock[j].id);
                }

                // Regardless of field-of-view, append j to i buffer list
                // if it is within the buffer region
                if (dist_ij < BUFFER_ALERT)
                    flock[i].buffer_list.push_back(flock[j].id);
            }
        }
}

void UpdateNeighboursFromBuffer(std::vector<Boid>& flock, std::vector<double>& dists)
{
    /**
        For each boid in the flock, evaluate the distance to other 
        boids within its buffer region and check if it's within close
        or near regions.
        If a boid is close (smallest radius category) then it is not
        included in the near list and vice versa - hence close boids
        only exert separation force, and near boids cohere and align.
        Buffer lists does not change.

        OUTPUTS
        Neighbour lists are stored within the Boid object as a class
        property.
        
        NOTES
        VISION_FOV, CLOSE_ALERT, NEAR_ALERT, defined in Settings.h.
        Must operate over all NxB boids where B is the number of boids
        in the buffer list of boid i.
    */

    int num_boids = flock.size();

    // Operate (parallel i) over all NxB boids - non-symmetric problem
    #pragma omp parallel for
    for (int i = 0; i < num_boids; i++)
        for (int id : flock[i].buffer_list)
        {
            // Read distance from ith to jth boid
            double dist_iid = dists[(i *num_boids) + id];

            // Calculate the angle from i velocity to id position
            float ang_iid = FovAngle(flock[i], flock[id]);

            // If id is within i field-of-view (i can see id) then 
            // append id to corresponding neighbour list of i
            if (ang_iid < VISION_FOV)
            {
                if (dist_iid < CLOSE_ALERT)
                    flock[i].close_list.push_back(flock[id].id);
                else if (dist_iid < NEAR_ALERT)
                    flock[i].near_list.push_back(flock[id].id);
            }
        }
}