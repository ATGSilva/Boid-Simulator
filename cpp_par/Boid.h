// -*- adamgillard-cpp -*-
#pragma once
#include "Vec3D.h"
#include <vector>
#include <memory>

class Boid
{
public:
    Boid();
    Boid(Vec3D, Vec3D, int);

    Vec3D pos, vel, acc;
    Vec3D cohere_force, sep_force, align_force, wind_force, wall_force;
    float mass;
    int id;

    // Neighbour values
    Vec3D sum_pos_cohere, sum_pos_sep, sum_vel_align;
    float sum_nmass, sum_cmass;
    int num_near, num_close;
};

void UpdatePos(Boid&);
void Reset(Boid&);

// Flock functions
std::vector<Boid> GenFlock(int, double, double, double, double);
