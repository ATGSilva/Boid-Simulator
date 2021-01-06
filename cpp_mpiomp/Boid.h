// -*- adamgillard-cpp -*-
#pragma once
#include "Vec3D.h"
#include <vector>
#include <memory>

struct Boid
{
    Boid();
    Boid(Vec3D, Vec3D, int);
    Boid(Vec3D, Vec3D, float, int);

    Vec3D pos;
    Vec3D vel;
    float mass;
    int id;
}__attribute__((__packed__));

struct Forces
{
    Forces();
    Vec3D cohere_force;
    Vec3D sep_force;
    Vec3D align_force;
    Vec3D wind_force;
    Vec3D wall_force;
};

void UpdatePos(Boid&, Forces&);

// Flock functions
std::vector<Boid> GenFlock(int, double, double, double, double);

