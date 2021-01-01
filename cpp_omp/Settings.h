// -*- adamgillard-cpp -*-
#pragma once
#include <cmath>

const short dims = 3;
const float duration = 200;
const float dt = 0.2;

const int buffer_for = 5;

const double wall_ubound = 1000;
const double wall_lbound = 0;

const double pos_ulim = wall_ubound;
const double pos_llim = wall_lbound;
const double vel_ulim = 15;
const double vel_llim = -10;

const float MAX_FORCE = 10;
const float MAX_ACC = 18;
const float MAX_VEL = 50;
const float WALL_FORCE = MAX_FORCE;
const float MAX_WIND = 3;

const float C_STR = 0.05;
const float S_STR = 0.25;
const float A_STR = 0.35;

const double buffer_alert = 300.0; // Defines buffer distance
const double near_alert = 250.0; // Defines near distance
const double close_alert = 10.0; // Defines close distance
const double vision_ang = (2.0f / 3.0f) * M_PI; // Defines angle over which a boid can see

