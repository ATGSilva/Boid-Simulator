/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Forces.h

FUNCTION SIGNATURE - RETURN TYPE
    CohereForce(const vector<Boid>&, Boid&, std::vector<int>&) - Vec3D
    SepForce(const vector<Boid>&, Boid&, std::vector<int>&) - Vec3D
    AlignForce(const vector<Boid>&, Boid&, std::vector<int>&) - Vec3D
    WallForce(Boid&) - Vec3D
    RandWindForce() - Vec3D
    WindEvo(Vec3D) - Vec3D
*/

#pragma once
#include <vector>
#include "Boid.h"
#include "Vec3D.h"

Vec3D CohereForce(std::vector<Boid>&, Boid&, std::vector<int>&);
Vec3D SepForce(std::vector<Boid>&, Boid&, std::vector<int>&);
Vec3D AlignForce(std::vector<Boid>&, Boid&, std::vector<int>&);
Vec3D WallForce(Boid&);
Vec3D RandWindForce();
Vec3D WindEvo(Vec3D);