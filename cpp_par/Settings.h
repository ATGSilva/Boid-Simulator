// -*- adamgillard-cpp -*-
#pragma once

const short dims = 3;
const float duration = 400;
const float dt = 0.5;

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
const float MAX_WIND = 4;

const float C_STR = 0.025;
const float S_STR = 0.35;
const float A_STR = 0.25;

