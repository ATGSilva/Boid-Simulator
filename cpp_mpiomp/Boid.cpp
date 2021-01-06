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
    mass = 1.f;
    id = ident;
}

Boid::Boid(Vec3D p, Vec3D v, float m, int ident)
{
    pos = p;
    vel = v;
    mass = m;
    id = ident;
}

Forces::Forces()
{
    Vec3D cohere_force;
    Vec3D sep_force;
    Vec3D align_force;
    Vec3D wind_force;
    Vec3D wall_force;
}

void UpdatePos(Boid& boid, Forces& force_list)
{
    Vec3D total_force;
    total_force = force_list.cohere_force + force_list.sep_force + force_list.align_force + force_list.wind_force;
    total_force = total_force.LimVec(MAX_FORCE);
    total_force = total_force + force_list.wall_force;
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
        rand_pos.z = 0.;
        Vec3D rand_vel = RandVec(vel_ulim, vel_llim);
        flock.emplace_back(Boid(rand_pos, rand_vel, i));
    }
    return flock;
}