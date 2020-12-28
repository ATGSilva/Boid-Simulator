// -*- adamgillard-cpp -*-
#pragma once


const int num_boids = 800;
const short dims = 3;
const float duration = 400;
const float dt = 0.5;

const double pos_ulim = 700;
const double pos_llim = 0;
const double vel_ulim = 15;
const double vel_llim = -10;

const float MAX_FORCE = 10;
const float MAX_ACC = 25;
const float MAX_VEL = 100;
const float WALL_FORCE = MAX_FORCE;

const float C_STR = 0.05;
const float S_STR = 0.35;
const float A_STR = 0.25;

const double wall_ubound = 800;
const double wall_lbound = 0;