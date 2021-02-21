/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Boid.cpp

Dependancy file to be used for building BoidSimOMP.cpp.
Contains boid struct and modification functions.

FUNCTION SIGNATURE - RETURN TYPE
    UpdatePos(Boid&) - void
    GenFlock(int, int, int, int, int) - std::vector<Boid>
*/

#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include <iostream>
#include <random>
#include <vector>
#include <memory>

// Boid generators ----------------------------------------------------
Boid::Boid() {}

Boid::Boid(Vec3D p, Vec3D v, int ident) 
{
    // Core properties only
    pos = p;
    vel = v;
    mass = 1.0f;
    id = ident;
}

Boid::Boid(Vec3D p, Vec3D v, float m, int ident)
{
    pos = p;
    vel = v;
    mass = m;
    id = ident;
}

// Force struct generator ---------------------------------------------
Forces::Forces()
{
    Vec3D cohere_force;
    Vec3D sep_force;
    Vec3D align_force;
    Vec3D wind_force;
    Vec3D wall_force;
}


void UpdatePos(Boid& boid, Forces& force_list, float Time)
{
    /**
        Updates the position of a boid using its forces. These are
        stored as struct properties so only input needed is the boid.
        Uses the Leapfrog method as position integrator. 
    */

    // Calculate total force as sum of components, then limit to max
    Vec3D total_force;
    total_force = force_list.cohere_force + force_list.sep_force + force_list.align_force + force_list.wind_force;
    total_force = total_force.LimVec(MAX_FORCE);
    total_force = total_force + force_list.wall_force;
    total_force = total_force.LimVec(MAX_FORCE);

    // Calculate the updates to acceleration, velocity and position
    Vec3D acc = ((1/boid.mass) * total_force).LimVec(MAX_ACC);
    // On initital update, offset velocity by half a timestep
    if (Time == 0)
    {
        boid.vel = (boid.vel + (DT/2)* acc).LimVec(MAX_VEL);
    } 
    // On every other timestep, update with offset carrying over
    else 
    {
        boid.vel = (boid.vel + (DT * acc)).LimVec(MAX_VEL);
    }
    boid.pos = boid.pos + (DT * boid.vel);
}

std::vector<Boid> GenFlock(int num_boids, int POS_ULIM, int POS_LLIM, int VEL_ULIM, int VEL_LLIM)
{
    /**
        Generates a standard library vector of length num_boids boids.
        Each boid is generated with random position and velocity.
        Z position set to 0 (taking off from the ground)
    */

    std::vector<Boid> flock;
    flock.reserve(num_boids);
    for (int i = 0; i < num_boids; i++) 
    {
        Vec3D rand_pos = RandVec(POS_ULIM, POS_LLIM);
        rand_pos.z = 0.0f;
        Vec3D rand_vel = RandVec(VEL_ULIM, VEL_LLIM);
        flock.emplace_back(Boid(rand_pos, rand_vel, i));
    }
    return flock;
}