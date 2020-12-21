#pragma once

const short MAX_VEL = 25;
const short MAX_ACC = 2;
const short MAX_FORCE = 10;
const short WALL_FORCE = MAX_FORCE;

const float C_STR = 0.4;
const float S_STR = 0.8;
const float A_STR = 0.5;

const int num_boids = 500;
const short dims = 3;
const float duration = 480;
const float dt = 0.1;

const double pos_ulim = 500;
const double pos_llim = 0;
const double vel_ulim = 10;
const double vel_llim = -vel_ulim;