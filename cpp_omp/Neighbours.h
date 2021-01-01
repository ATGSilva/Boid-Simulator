#pragma once
#include <vector>
#include <memory>
#include "Boid.h"

std::vector<double> FindDists(std::vector<Boid>&);
void FindNeighbours(std::vector<Boid>&, std::vector<double>&);
void UpdateNeighboursFromBuffer(std::vector<Boid>&, std::vector<double>&);