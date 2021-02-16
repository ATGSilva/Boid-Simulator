/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Boid.h

Dependancy file to be used in BoidSimOMP.cpp.
Contains boid struct and modification functions.

FUNCTION SIGNATURE - RETURN TYPE
    Reset(Boid&, int) - void
    UpdatePos(Boid&) - void
    GenFlock(int, int, int, int, int) - std::vector<Boid>
*/

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

void Reset(Boid&, int);
void UpdatePos(Boid&, float);

// Flock functions
std::vector<Boid> GenFlock(int, int, int, int, int);
