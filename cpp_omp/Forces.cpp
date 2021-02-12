/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Forces.cpp

Dependancy file to be used for building BoidSimOMP.cpp.
Contains boid force calculation functions.

FUNCTION SIGNATURE - RETURN TYPE
    CohereForce(const vector<Boid>&, Boid&) - void
    SepForce(const vector<Boid>&, Boid&) - void
    AlignForce(const vector<Boid>&, Boid&) - void
    NormaliseForce(Vec3D) - void
    WallForce(Boid&) - void
    RandWindForce - Vec3D
    WindEvo(Vec3D) - Vec3D
*/

#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <memory>
#include <omp.h>
#include "Boid.h"
#include "Forces.h"
#include "Vec3D.h"
#include "Settings.h"


// Core behavioural forces --------------------------------------------
void CohereForce(const std::vector<Boid>& flock, Boid& boid)
{
    /**
        Calculates the coherence pseudo-force acting on a single boid 
        from all "near" boids.
        If no boids are "near" then returns a zero vector.

        EQUATION
        Force = [(1/sum_neighbour_masses) * COM_neighbours] - boid.position

        NOTES
        Force is arbitrarily scaled in magnitude to exaggerate behaviour,
        this can be done because we are using a pseudo-force.
    */

    int num_near = boid.near_list.size();
    // If neighbours exist, calculate force, else return zero force vector
    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bpos = boid.pos;
        const Vec3D bvel = boid.vel;

        // For every id in near boid list, calculate mass sum and com 
        for (int id : boid.near_list)
        {
            m_sum += flock[id].mass;
            pos_sum = pos_sum + (flock[id].mass * flock[id].pos);
        }

        Vec3D cohere_force = ((1/m_sum) * pos_sum) - bpos;  // Calculate Force
        NormaliseForce(cohere_force);                       // Normalise Force (see function)
        boid.cohere_force = cohere_force * C_STR;           // Write scaled force to boid properties
    }
    else boid.cohere_force = Vec3D(); 
}

void SepForce(const std::vector<Boid>& flock, Boid& boid)
{
    /**
        Calculates the separation pseudo-force acting on a single boid 
        from all "close" boids.
        If no boids are "close" then returns a zero vector.

        EQUATION
        Force = (1/sum_neighbour_masses) * Centre of Displacements

        NOTES
        Force is arbitrarily scaled in magnitude to exaggerate behaviour,
        this can be done because we are using a pseudo-force.
    */

    int num_close = boid.close_list.size();
    // If neighbours exist, calculate force, else return zero force vector
    if (num_close != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        // For every id in close boid list, calculate mass sum and centre of displacements 
        for (int id : boid.close_list)
        {
            m_sum += flock[id].mass;
            Vec3D disp_idb = boid.pos - flock[id].pos;
            pos_sum = pos_sum + (flock[id].mass * disp_idb);
        }

        Vec3D sep_force = (1/m_sum) * pos_sum;  // Calculate Force
        NormaliseForce(sep_force);              // Normalise force
        boid.sep_force = sep_force * S_STR;     // Write scaled force to boid properties
    }
    else boid.sep_force = Vec3D();
}

void AlignForce(const std::vector<Boid>& flock, Boid& boid)
{
    /**
        Calculates the alignment pseudo-force acting on a single boid 
        from all "near" boids.
        If no boids are "near" then returns a zero vector.

        EQUATION
        Force = [(1/sum(neighbour_masses)) * sum(neighbour_velcoities)] - boid.velocity

        NOTES
        Force is arbitrarily scaled in magnitude to exaggerate behaviour,
        this can be done because we are using a pseudo-force.
    */

    int num_near = boid.near_list.size();

    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D vel_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        // For every id in close boid list, calculate mass sum and velocity sum
        for (int id : boid.near_list)
        {
            m_sum += flock[id].mass;
            vel_sum = vel_sum + (flock[id].mass * flock[id].vel);
        }

        Vec3D align_force = (1/m_sum) * vel_sum;            // Calculate Force
        NormaliseForce(align_force);                        // Normalise Force
        boid.align_force = (align_force - bvel) * A_STR;    // Write scaled force to boid properties
    }
    else boid.align_force = Vec3D();
}

void NormaliseForce(Vec3D& force)
{
    /**
        Normalises a given force to the maximum velocity value.
    */

    double force_mag = sqrt(pow(force.x, 2) + pow(force.y, 2) + pow(force.z, 2));
    if (force_mag != 0)
        force = force * (MAX_VEL / force_mag);
}

// Bounding wall force ------------------------------------------------
void WallForce(Boid& boid)
{
    /**
        Calculates the bounding force acting on a single boid.
        The force scales with the distance of the boid from the
        bounding wall position in each cartesian axis.
        Bounding box values are provided in Settings.h.
    */

    const double x = boid.pos.x;
    const double y = boid.pos.y;
    const double z = boid.pos.z;

    int x_ubound = (x > WALL_UBOUND);
    int y_ubound = (y > WALL_UBOUND);
    int z_ubound = (z > WALL_UBOUND);

    int x_lbound = (x < WALL_LBOUND);
    int y_lbound = (y < WALL_LBOUND);
    int z_lbound = (z < WALL_LBOUND);

    Vec3D which_pos_wall = Vec3D(x_ubound, y_ubound, z_ubound);
    Vec3D which_neg_wall = Vec3D(x_lbound, y_lbound, z_lbound);
    Vec3D u_bounds = Vec3D(WALL_UBOUND, WALL_UBOUND, WALL_UBOUND);
    Vec3D l_bounds = Vec3D(WALL_LBOUND, WALL_LBOUND, WALL_LBOUND);
    Vec3D wall_force = Vec3D();

    if (x_ubound || y_ubound || z_ubound)
        wall_force = wall_force + (-1) * WALL_FORCE * (which_pos_wall * (boid.pos - u_bounds));

    if (x_lbound || y_lbound || z_lbound)
        wall_force = wall_force +  WALL_FORCE * (which_neg_wall * (l_bounds - boid.pos));

    boid.wall_force = wall_force.LimVec(WALL_FORCE);
}

// Wind force functions -----------------------------------------------
Vec3D RandWindForce()
{
    /**
        Generates a random wind force vector where values are
        between a supplied (Settings.h) upper and lower bound.
    */

    double uwind_force = MAX_WIND;
    double lwind_force = -MAX_WIND;
    return RandVec(lwind_force, uwind_force);
}

Vec3D WindEvo(Vec3D wind_force)
{
    /**
        Randomly evolves a wind force.
        The magnitude of the evolution is limited to between the +- 
        value of the timestep multiplied by a random factor between 0.1
        and 4, creating a more realistic wind evolution.
        Wind force is limited to a maximum value defined in Settings.h.

    */

    double uwind_change = DT;
    double lwind_change = -DT;
    double rand_factor = RandVal(0.1, 4);
    Vec3D wind_change = rand_factor * RandVec(lwind_change, uwind_change);

    return (wind_force + wind_change).LimVec(MAX_WIND);
}