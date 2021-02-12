#pragma once
#include <chrono>
#include <fstream>
#include <sstream>
#include <omp.h>
#include "Boid.h"

void BeginTimingSession(std::ofstream& timingfile, const int& num_boids, const int& threads)
{
    std::ostringstream tname;
    std::string tfilename;

    tname << "data/" << "B" << num_boids << "Thr" << threads << ".csv";
    tfilename = tname.str();
    timingfile.open(tfilename);
    timingfile << "Frame,Function,Duration\n";
}

void BeginResultSession(std::ofstream& results)
{
    results.open("pos_dat.csv");
    results << "Time,Frame,ID,Mass,X,Y,Z,VelX,VelY,VelZ\n";
}

void WriteResults(std::ofstream& results, float time, int iter, Boid& boid)
{
    results << time << "," << iter << "," << boid.id << "," << boid.mass << "," << boid.pos << "," << boid.vel << "\n";
}

void EndWriteSession(std::ofstream& file)
{
    file.close();
}

void WriteTiming(std::ofstream& timingfile, const char* func_name, int iters, std::chrono::microseconds duration)
{
    timingfile << iters << "," << func_name << "," << duration.count()/1e6 << "\n";
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