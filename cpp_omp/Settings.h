/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Settings.h

Contains constants to be used for files:
    BoidSimOMP.cpp
    BoidSimOMPMPI.cpp
    Neighbours.cpp
    Forces.cpp
    Boid.cpp
    Vec3D.cpp
*/
#pragma once
// Include core header files
#include <cmath>

// Timestep Settings ---------------------
const float DURATION = 400.0f;                      // Total duration (seconds)
const float DT = 0.5f;                              // Timestep (seconds)
// Buffer Settings -----------------------
const bool BUFFER = false;                          // Turn buffer on or off
const int BUFFER_FOR = 1;                           // Number of iterations to buffer over
// Bound Limit Settings ------------------
const int WALL_UBOUND = 1000;                       // Upper wall boundary in all dimensions
const int WALL_LBOUND = 0;                          // Lower wall boundary in all dimensions
const int POS_ULIM = WALL_UBOUND;                   // Upper initial position limit
const int POS_LLIM = WALL_LBOUND;                   // Lower initial position limit
const int VEL_ULIM = 15;                            // Upper initial velocity limit
const int VEL_LLIM = -10;                           // Lower initial velocity limit
// Force Bound Settings ------------------
const float MAX_FORCE = 10.0f;                      // Maximum force on boid
const float MAX_ACC = 18.0f;                        // Maximum  boid acceleration
const float MAX_VEL = 60.0f;                        // Maximum boid velocity
const float WALL_FORCE = MAX_FORCE;                 // Wall force strength
const float MAX_WIND = 1.5f;                        // Maximum wind force
// Force Strength Multipliers ------------
const float C_STR = 0.06f;                          // Coherence force strength modifier
const float S_STR = 0.25f;                          // Separation force strength modifier
const float A_STR = 0.35f;                          // Alignment force strength modifier
// Neighbour Finding Settings ------------
const int BUFFER_ALERT = 300;                       // Defines buffer distance
const int NEAR_ALERT = 250;                         // Defines near distance
const int CLOSE_ALERT = 10;                         // Defines close distance
const float VISION_FOV = (2.0f / 3.0f) * M_PI;      // Defines angle over which a boid can see
// File Output Control -------------------
//const int TIMING = 1;                           // If true, readout timing data to file
//const int BOID_READOUT = 1;                     // If true, readout boid property data to file
