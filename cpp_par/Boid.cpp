// -*- adamgillard-cpp -*-
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include <iostream>
#include <random>
#include <vector>
#include <memory>


Boid::Boid()
{
    pos = Vec3D();
    vel = Vec3D();
    mass = 1;
    cohere_force = Vec3D();
    sep_force = Vec3D();
    align_force = Vec3D();
    wind_force = Vec3D();
    wall_force = Vec3D();
    sum_pos_cohere = Vec3D();
    sum_pos_sep = Vec3D();
    sum_vel_align = Vec3D();
    sum_nmass = 0;
    sum_cmass = 0;
    num_near = 0;
    num_close = 0;
}

Boid::Boid(Vec3D p, Vec3D v, int ident) 
{
    pos = p;
    vel = v;
    mass = 1;
    id = ident;
    cohere_force = Vec3D();
    sep_force = Vec3D();
    align_force = Vec3D();
    wind_force = Vec3D();
    wall_force = Vec3D();
    sum_pos_cohere = Vec3D();
    sum_pos_sep = Vec3D();
    sum_vel_align = Vec3D();
    sum_nmass = 0;
    sum_cmass = 0;
    num_near = 0;
    num_close = 0;
}

void Reset(Boid& boid)
{   
    boid.cohere_force = Vec3D();
    boid.sep_force = Vec3D();
    boid.align_force = Vec3D();
    boid.wind_force = Vec3D();
    boid.wall_force = Vec3D();
    boid.sum_pos_cohere = Vec3D();
    boid.sum_pos_sep = Vec3D();
    boid.sum_vel_align = Vec3D();
    boid.sum_nmass = 0;
    boid.sum_cmass = 0;
    boid.num_near = 0;
    boid.num_close = 0;
}

void UpdatePos(Boid& boid)
{
    Vec3D total_force;
    total_force = boid.cohere_force + boid.sep_force + boid.align_force + boid.wind_force;
    total_force = total_force.LimVec(MAX_FORCE);
    total_force = total_force + boid.wall_force;
    total_force = total_force.LimVec(MAX_FORCE);

    Vec3D acc = ((1/boid.mass) * total_force).LimVec(MAX_ACC);
    boid.vel = (boid.vel + (dt * acc)).LimVec(MAX_VEL);
    boid.pos = boid.pos + (dt * boid.vel);
}

std::vector<Boid> GenFlock(int num_boids, double pos_ulim, double pos_llim, double vel_ulim, double vel_llim)
{
    std::vector<Boid> flock;
    flock.reserve(num_boids);
    for (int i = 0; i < num_boids; i++) 
    {
        Vec3D rand_pos = RandVec(pos_ulim, pos_llim);
        Vec3D rand_vel = RandVec(vel_ulim, vel_llim);
        flock.emplace_back(rand_pos, rand_vel, i);
    }
    return flock;
}