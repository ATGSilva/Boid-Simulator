/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Settings.h

Contains constants to be used for files:
    BoidSimOMP.cpp
    Neighbours.cpp
    Forces.cpp
    Boid.cpp
    Vec3D.cpp
*/
#pragma once
// Include core header files
#include <cmath>

// Timestep Settings ---------------------
const float DURATION = 400.0f;                      // Total duration (seconds) (400)
const float DT = 0.5f;                              // Timestep (seconds) (0.5)
// Buffer Settings -----------------------
const bool BUFFER = false;                          // Turn buffer on or off
const int BUFFER_FOR = 1;                           // Number of iterations to buffer over (1 = no buffering)
// Bound Limit Settings ------------------
const int WALL_UBOUND = 1000;                       // Upper wall boundary in all dimensions (1000)
const int WALL_LBOUND = 0;                          // Lower wall boundary in all dimensions (0)
const int POS_ULIM = WALL_UBOUND;                   // Upper initial position limit
const int POS_LLIM = WALL_LBOUND;                   // Lower initial position limit
const int VEL_ULIM = 45;                            // Upper initial velocity limit (45)
const int VEL_LLIM = -20;                           // Lower initial velocity limit (-20)
// Force Bound Settings ------------------
const float MAX_FORCE = 10.0f;                      // Maximum force on boid (10)
const float MAX_ACC = 18.0f;                        // Maximum  boid acceleration (18)
const float MAX_VEL = 50.0f;                        // Maximum boid velocity (50)
const float WALL_FORCE = MAX_FORCE;                 // Wall force strength (MAX_FORCE)
const float MAX_WIND = 1.5f;                        // Maximum wind force (1.5)
// Force Strength Multipliers ------------
const float C_STR = 0.06f;                          // Coherence force strength modifier (0.06)
const float S_STR = 0.25f;                          // Separation force strength modifier (0.25)
const float A_STR = 0.35f;                          // Alignment force strength modifier (0.35)
// Neighbour Finding Settings ------------
const int BUFFER_ALERT = 300;                       // Defines buffer distance (300)
const int NEAR_ALERT = 200;                         // Defines near distance (200)
const int CLOSE_ALERT = 10;                         // Defines close distance (10)
const float VISION_FOV = (2.0f / 3.0f) * M_PI;      // Defines angle over which a boid can see (2/3 pi) (pi = can see all around)
// File Output Control -------------------
#define TIMING 1                                    // If 1, readout timing data to file, else 0
#define BOID_READOUT 1                              // If 1, readout boid property data to file, else 0
#define PROGRESS 1                                  // If 1, show progress bar updates, else 0
#define VISION 1
