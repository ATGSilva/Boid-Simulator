#pragma once

const int num_boids = 1000;
const short dims = 3;
const float duration = 400;
const float dt = 0.1;

const double pos_ulim = 300;
const double pos_llim = 0;
const double vel_ulim = 15;
const double vel_llim = -10;

const float MAX_VEL = 25;
const float MAX_ACC = 5;
const float MAX_FORCE = 10;
const float WALL_FORCE = 5;

const float C_STR = 0.25;
const float S_STR = 0.8;
const float A_STR = 0.5;

const double wall_ubound = 1500;
const double wall_lbound = 0;