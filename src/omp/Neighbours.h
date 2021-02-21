/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Neighbours.cpp

Dependancy file to be used for building BoidSimOMP.cpp

FUNCTION SIGNATURE - RETURN TYPE
    FindDists(vector<Boid>&) - vector<double>
    FindNeighbours(vector<Boid>&, vector<double>&) - void
    UpdateNeighboursFromBuffer(vector<Boid>&, vector<double>&) - void
*/

#pragma once
#include <vector>
#include <memory>
#include "Boid.h"

std::vector<double> FindDists(std::vector<Boid>&);
void FindNeighbours(std::vector<Boid>&, std::vector<double>&);
void UpdateNeighboursFromBuffer(std::vector<Boid>&, std::vector<double>&);