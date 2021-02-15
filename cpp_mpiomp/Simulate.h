/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Simulate.h

Dependancy file to be used for building BoidSimMPIOMP.cpp
Contains boid neighbour finding routines and simulation wrapper.

FUNCTION SIGNATURE - RETURN TYPE
    FindDists(vector<Boid>&, int, int, int) - vector<double>
    FindNeighbours(vector<Boid>&, int, vector<double>&) - void
    Simulate(std::vector<Boid>&, int, int, std::vector<double>&, Vec3D, float)
*/

#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "Boid.h"

std::vector<double> FindDists(std::vector<Boid>&, int, int, int);
std::vector<Boid> Simulate(std::vector<Boid>&, int, int, std::vector<double>&, Vec3D);


