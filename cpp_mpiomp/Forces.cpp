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


Vec3D CohereForce(std::vector<Boid>& flock, Boid& boid, std::vector<int>& near_list)
{
    Vec3D cohere_force;
    int num_near = near_list.size();
    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bpos = boid.pos;
        const Vec3D bvel = boid.vel;

        for (int id : near_list)
        {
            m_sum += flock[id].mass;
            pos_sum = pos_sum + (flock[id].mass * flock[id].pos);
        }

        Vec3D target = ((1/m_sum) * pos_sum) - bpos;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        cohere_force = target * C_STR;
    }
    else cohere_force = Vec3D();

    return cohere_force;
}

Vec3D SepForce(std::vector<Boid>& flock, Boid& boid, std::vector<int>& close_list)
{
    Vec3D sep_force;
    int num_close = close_list.size();
    if (num_close != 0)
    {
        float m_sum = 0;
        Vec3D pos_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        for (int id : close_list)
        {
            Vec3D disp_idb = boid.pos - flock[id].pos;
            m_sum += flock[id].mass;
            pos_sum = pos_sum + (flock[id].mass * disp_idb);
        }

        Vec3D target = (1/m_sum) * pos_sum;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        sep_force = target * S_STR;
    }
    else sep_force = Vec3D();

    return sep_force;
}

Vec3D AlignForce(std::vector<Boid>& flock, Boid& boid, std::vector<int>& near_list)
{
    Vec3D align_force;
    int num_near = near_list.size();
    if (num_near != 0)
    {
        float m_sum = 0;
        Vec3D vel_sum = Vec3D();
        const Vec3D bvel = boid.vel;

        for (int id : near_list)
        {
            m_sum += flock[id].mass;
            vel_sum = vel_sum + (flock[id].mass * flock[id].vel);
        }

        Vec3D target = (1/m_sum) * vel_sum;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        align_force = (target - bvel) * A_STR; 
    }
    else align_force = Vec3D();

    return align_force;
}

Vec3D WallForce(Boid& boid)
{
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

    return wall_force.LimVec(WALL_FORCE);
}

Vec3D RandWindForce()
{
    double uwind_force = MAX_WIND;
    double lwind_force = -MAX_WIND;
    return RandVec(lwind_force, uwind_force);
}

Vec3D WindEvo(Vec3D wind_force)
{
    double uwind_change = DT;
    double lwind_change = -DT;
    double rand_factor = RandVal(0.1, 4);
    Vec3D wind_change = rand_factor * RandVec(lwind_change, uwind_change);

    return (wind_force + wind_change).LimVec(MAX_WIND);
}

