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
    std::vector<int> near_list, close_list, buffer_list;
};

void UpdatePos(Boid&);
void Reset(Boid&, int);

// Flock functions
std::vector<Boid> GenFlock(int, double, double, double, double);
