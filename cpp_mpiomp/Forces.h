// -*- adamgillard-cpp -*-
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