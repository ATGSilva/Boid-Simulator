#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "Boid.h"

void ProgBar(float&);
std::vector<double> FindDists(std::vector<Boid>&, int, int, int);
std::vector<Boid> Simulate(std::vector<Boid>&, int, int, std::vector<double>&, Vec3D, float);


