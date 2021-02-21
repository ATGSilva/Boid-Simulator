/**
-*- adamgillard-cpp -*- Advanced Computational Physics -*-

Timing.h

Dependancy file to be used in BoidSimOMP.cpp.
Contains functions to write data to file, including timings and boid
properties.

FUNCTION SIGNATURE - RETURN TYPE
    BeginTimingSession(std::ofstream&, const int&, const int&) - void
    BeginResultSession(std::ofstream&) - void
    WriteResults(std::ofstream&, float, int, Boid&) - void
    WriteTiming(std::ofstream&, const char*, int, std::chrono::microseconds) - void
    EndWriteSession(std::ofstream&) - void
*/

#pragma once
#include <chrono>
#include <fstream>
#include <sstream>
#include <omp.h>

void BeginTimingSession(std::ofstream& timingfile, const int& num_boids, const int& threads)
{
    std::ostringstream tname;
    std::string tfilename;

    tname << "timing/" << "B" << num_boids << "Thr" << threads << ".csv";
    tfilename = tname.str();
    timingfile.open(tfilename);
    timingfile << "Frame,Function,Duration\n";
}

void BeginResultSession(std::ofstream& results, const int& num_boids, const int& threads)
{
    std::ostringstream rname;
    std::string rfilename;

    rname << "boid_data/" << "B" << num_boids << "Thr" << threads << ".csv";
    rfilename = rname.str();
    results.open(rfilename);
    results << "Time,Frame,ID,Mass,X,Y,Z,VelX,VelY,VelZ\n";
}

void WriteResults(std::ofstream& results, float time, int iter, Boid& boid)
{
    results << time << "," << iter << "," << boid.id << "," << boid.mass << "," << boid.pos << "," << boid.vel << "\n";
}

void WriteTiming(std::ofstream& timingfile, const char* func_name, int iters, std::chrono::microseconds duration)
{
    timingfile << iters << "," << func_name << "," << duration.count()/1e6 << "\n";
}

void EndWriteSession(std::ofstream& file)
{
    file.close();
}

class Timer
{
public:
    Timer(const char* name, std::ofstream& file, int iter) : m_Name(name), m_Stopped(false), timingfile(file)
    {
        tstart = std::chrono::high_resolution_clock::now();
        iters = iter;
    };

    ~Timer()
    {
        if (!m_Stopped)
            Stop();
    }

    void Stop()
    {
        auto tend = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);

        #pragma omp critical
        WriteTiming(timingfile, m_Name, iters, duration);

        m_Stopped = true;
    }
private:
    const char* m_Name;
    std::chrono::time_point<std::chrono::high_resolution_clock> tstart;
    std::ofstream& timingfile;
    int iters;
    bool m_Stopped;
};

void WriteSettings(int num_boids, int threads)
{
    std::ofstream settings;
    settings.open("boid_data/settings.cfg");

    settings << "[Program]\n";
    settings << "type = OpenMP\n";
    settings << "number_boids = " << num_boids << "\n";
    settings << "number_procs = " << threads << "\n";
    settings << "\n[Time]\n";
    settings << "duration = " << DURATION << "\n";
    settings << "timestep = " << DT << "\n";
    settings << "\n[Buffer]\n";
    settings << "buffer = " << BUFFER << "\n";
    settings << "buffer_for = " << BUFFER_FOR << "\n";
    settings << "\n[Bounds]\n";
    settings << "wall_bound_upper = " << WALL_UBOUND << "\n";
	settings << "wall_bound_lower = " << WALL_LBOUND << "\n";
    settings << "random_position_bounds = [" << POS_LLIM << ", " << POS_ULIM << "]\n";
    settings << "random_velocity_bounds = [" << VEL_LLIM << ", " << VEL_ULIM << "]\n";
    settings << "\n[Kinematic Limits]\n";
    settings << "max_force = " << MAX_FORCE << "\n";
    settings << "max_acceleration = " << MAX_ACC << "\n";
    settings << "max_velocity = " << MAX_VEL << "\n";
    settings << "wall_force = " << WALL_FORCE << "\n";
    settings << "max_wind_force = " << MAX_WIND << "\n";
    settings << "\n[Force Multipliers]\n";
    settings << "coherence_multiplier = " << C_STR << "\n";
    settings << "separation_multipler = " << S_STR << "\n";
    settings << "allignment_multipler = " << A_STR << "\n";
    settings << "\n[Neighbour Finding]\n";
    settings << "buffer_radius = " << BUFFER_ALERT << "\n";
    settings << "near_radius = " << NEAR_ALERT << "\n";
    settings << "too_close_radius = " << CLOSE_ALERT << "\n";
    settings << "field_of_view = " << VISION_FOV << "\n";
    settings << "\n[Control]\n";
    settings << "timing_readout = " << TIMING << "\n";
    settings << "boid_readout = " << BOID_READOUT << "\n";
    settings << "progress_bar_onoff = " << PROGRESS << "\n";
    settings << "vision_onoff = " << VISION << "\n";

    settings.close();
}