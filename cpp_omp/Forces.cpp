// -*- adamgillard-cpp -*-
#include "Boid.h"
#include "Forces.h"
#include "Vec3D.h"
#include "Settings.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <memory>


void CohereForce(const std::vector<Boid>& flock, Boid& boid)
{
    int num_near = boid.near_list.size();
    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bpos = boid.pos;
        const Vec3D bvel = boid.vel;

        for (int id : boid.near_list)
        {
            m_sum += flock[id].mass;
            pos_sum = pos_sum + ((1/flock[id].mass) * flock[id].pos);
        }

        Vec3D target = ((1/m_sum) * pos_sum) - bpos;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        boid.cohere_force = (target - bvel) * C_STR;
    }
    else boid.cohere_force = Vec3D(); 
}

void SepForce(const std::vector<Boid>& flock, Boid& boid)
{
    int num_close = boid.close_list.size();
    if (num_close != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        for (int id : boid.close_list)
        {
            Vec3D disp_idb = boid.pos - flock[id].pos;
            m_sum += flock[id].mass;
            pos_sum = pos_sum + (flock[id].mass * disp_idb);
        }

        Vec3D target = (1/m_sum) * pos_sum;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        boid.sep_force = (target - bvel) * S_STR;
    }
    else boid.sep_force = Vec3D();
}

void AlignForce(const std::vector<Boid>& flock, Boid& boid)
{
    int num_near = boid.near_list.size();

    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D vel_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        for (int id : boid.near_list)
        {
            m_sum += flock[id].mass;
            vel_sum = vel_sum + (flock[id].mass * flock[id].vel);
        }

        Vec3D target = (1/m_sum) * vel_sum;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        boid.align_force = (target - bvel) * A_STR; 
    }
    else boid.align_force = Vec3D();
}

void WallForce(Boid& boid)
{
    const double x = boid.pos.x;
    const double y = boid.pos.y;
    const double z = boid.pos.z;

    int x_ubound = (x > wall_ubound);
    int y_ubound = (y > wall_ubound);
    int z_ubound = (z > wall_ubound);

    int x_lbound = (x < wall_lbound);
    int y_lbound = (y < wall_lbound);
    int z_lbound = (z < wall_lbound);

    Vec3D which_pos_wall = Vec3D(x_ubound, y_ubound, z_ubound);
    Vec3D which_neg_wall = Vec3D(x_lbound, y_lbound, z_lbound);
    Vec3D u_bounds = Vec3D(wall_ubound, wall_ubound, wall_ubound);
    Vec3D l_bounds = Vec3D(wall_lbound, wall_lbound, wall_lbound);
    Vec3D wall_force = Vec3D();

    if (x_ubound || y_ubound || z_ubound)
        wall_force = wall_force + (-1) * WALL_FORCE * (which_pos_wall * (boid.pos - u_bounds));

    if (x_lbound || y_lbound || z_lbound)
        wall_force = wall_force +  WALL_FORCE * (which_neg_wall * (l_bounds - boid.pos));

    boid.wall_force = wall_force.LimVec(WALL_FORCE);
}

Vec3D RandWindForce()
{
    double uwind_force = MAX_WIND;
    double lwind_force = -MAX_WIND;
    return RandVec(lwind_force, uwind_force);
}

Vec3D WindEvo(Vec3D wind_force)
{
    double uwind_change = dt;
    double lwind_change = -dt;
    double rand_factor = RandVal(0.1, 4);
    Vec3D wind_change = rand_factor * RandVec(lwind_change, uwind_change);

    return (wind_force + wind_change).LimVec(MAX_WIND);
}

