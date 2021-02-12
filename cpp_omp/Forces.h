/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Forces.h

FUNCTION SIGNATURE - RETURN TYPE
    CohereForce(const vector<Boid>&, Boid&) - void
    SepForce(const vector<Boid>&, Boid&) - void
    AlignForce(const vector<Boid>&, Boid&) - void
    WallForce(Boid&) - void
    RandWindForce - Vec3D
    WindEvo(Vec3D) - Vec3D
*/
#pragma once

void CohereForce(const std::vector<Boid>&, Boid&);
void SepForce(const std::vector<Boid>&, Boid&);
void AlignForce(const std::vector<Boid>&, Boid&);
void WallForce(Boid&);
Vec3D RandWindForce();
Vec3D WindEvo(Vec3D);