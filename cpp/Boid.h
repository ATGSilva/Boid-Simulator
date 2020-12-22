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
    Vec3D cohere_force, sep_force, align_force, wall_force;
    float mass;
    int id;

    // Neighbour values
    Vec3D sum_pos_cohere;
    Vec3D sum_pos_sep;
    Vec3D sum_vel_align;
    float sum_nmass;
    float sum_cmass;
    int num_near;
    int num_close;
};

void UpdatePos(Boid&);
void Reset(Boid&);

// Flock functions
std::vector<Boid> GenFlock(int, double, double, double, double);
std::vector<Boid> ResetBoids(std::vector<Boid>);
