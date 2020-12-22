#include "Boid.h"
#include "Forces.h"
#include "Vec3D.h"
#include "Settings.h"
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include <memory>


void CohereForce(Boid& boid)
{
    int num_near = boid.num_near;
    if (num_near != 0)
    {
        const float m_sum = boid.sum_nmass;
        const Vec3D pos_sum = boid.sum_pos_cohere;
        const Vec3D bpos = boid.pos;
        const Vec3D bvel = boid.vel;

        Vec3D target = ((1/m_sum) * pos_sum) - bpos;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        boid.cohere_force = (target - bvel) * C_STR;
    }
    else boid.cohere_force = Vec3D(); 
}

void SepForce(Boid& boid)
{
    int num_close = boid.num_close;
    if (num_close != 0)
    {
        const float m_sum = boid.sum_nmass;
        const Vec3D pos_sum = boid.sum_pos_sep;
        const Vec3D bvel = boid.vel;

        Vec3D target = (1/m_sum) * pos_sum;
        double target_mag = sqrt(pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2));
        if (target_mag != 0)
            target = target * (MAX_VEL / target_mag);

        boid.sep_force = (target - bvel) * S_STR;
    }
    else boid.sep_force = Vec3D();
}

void AlignForce(Boid& boid)
{
    int num_near = boid.num_near;

    if (num_near != 0)
    {
        const float m_sum = boid.sum_nmass;
        const Vec3D vel_sum = boid.sum_vel_align;
        const Vec3D bvel = boid.vel;

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
        wall_force = wall_force + WALL_FORCE * (which_neg_wall * (l_bounds - boid.pos));

    boid.wall_force = wall_force.LimVec(WALL_FORCE);
}