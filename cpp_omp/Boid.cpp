// -*- adamgillard-cpp -*-
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include <iostream>
#include <random>
#include <vector>
#include <memory>


Boid::Boid() {}

Boid::Boid(Vec3D p, Vec3D v, int ident) 
{
    pos = p;
    vel = v;
    mass = 1;
    id = ident;

    // Properties to reset on each timestep
    cohere_force = Vec3D();
    sep_force = Vec3D();
    align_force = Vec3D();
    wind_force = Vec3D();
    wall_force = Vec3D();

    near_list = std::vector<int>();
    close_list = std::vector<int>();
    buffer_list = std::vector<int>();
}

void Reset(Boid& boid, int iters)
{   
    boid.cohere_force = Vec3D();
    boid.sep_force = Vec3D();
    boid.align_force = Vec3D();
    boid.wind_force = Vec3D();
    boid.wall_force = Vec3D();
    boid.near_list.clear();
    boid.close_list.clear();
    if (iters % BUFFER_FOR != 0 || BUFFER_FOR == 1)
        boid.buffer_list.clear();
}

void UpdatePos(Boid& boid)
{
    Vec3D total_force;
    total_force = boid.cohere_force + boid.sep_force + boid.align_force + boid.wind_force;
    total_force = total_force.LimVec(MAX_FORCE);
    total_force = total_force + boid.wall_force;
    total_force = total_force.LimVec(MAX_FORCE);

    Vec3D acc = ((1/boid.mass) * total_force).LimVec(MAX_ACC);
    boid.vel = (boid.vel + (DT * acc)).LimVec(MAX_VEL);
    boid.pos = boid.pos + (DT * boid.vel);
}

std::vector<Boid> GenFlock(int num_boids, int POS_ULIM, int POS_LLIM, int VEL_ULIM, int VEL_LLIM)
{
    std::vector<Boid> flock;
    flock.reserve(num_boids);
    for (int i = 0; i < num_boids; i++) 
    {
        Vec3D rand_pos = RandVec(POS_ULIM, POS_LLIM);
        rand_pos.z = 0;
        Vec3D rand_vel = RandVec(VEL_ULIM, VEL_LLIM);
        flock.emplace_back(rand_pos, rand_vel, i);
    }
    return flock;
}