// -*- adamgillard-cpp -*-
#pragma once

void CohereForce(const std::vector<Boid>&, Boid&);
void SepForce(const std::vector<Boid>&, Boid&);
void AlignForce(const std::vector<Boid>&, Boid&);
void WallForce(Boid&);
Vec3D RandWindForce();
Vec3D WindEvo(Vec3D);