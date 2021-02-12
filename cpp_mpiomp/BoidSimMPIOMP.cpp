// -*- adamgillard-cpp -*-

// COMPILE: mpicxx -fopenmp BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim -O3
// RUN: mpirun -np <NUM_PROCS> ./BoidSim <THREADS_PER_PROC> <NUM_BOIDS>
#include "Boid.h"
#include "Vec3D.h"
#include "Forces.h"
#include "Settings.h"
#include "Simulate.h"
#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <memory>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <mpi.h>

int ArgChecks(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: ./BoidSimOMP <Integer NUM_THREADS> <Integer NUM_BOIDS>" << std::endl;
        return 1;
    }

    int threads = std::strtol(argv[1], nullptr, 0);
    if (threads < 1 || threads > omp_get_max_threads())
    {
        std::cerr << "<NUM_THREADS> must be between 1 and " << omp_get_max_threads() << std::endl;
        return 1;
    }
    int num_boids = std::strtol(argv[2], nullptr, 0);
    if (num_boids < 2)
    {
        std::cerr << "<NUM_BOIDS> must be a positive integer more than 1." << std::endl;
        return 1;
    }
    else return 0;
}

MPI_Datatype VecDataType()
{
    Vec3D vec = Vec3D(1., 2., 3.);
    MPI_Datatype VECTYPE;
    int vlens[3] = {1,1,1};
    MPI_Aint vdisps[3] = {0};
    MPI_Aint vst_addr, vaddr;
    MPI_Get_address(&vec.x, &vst_addr);
    MPI_Get_address(&vec.y, &vaddr);
    vdisps[1] = vaddr - vst_addr;
    MPI_Get_address(&vec.z, &vaddr);
    vdisps[2] = vaddr - vst_addr;
    MPI_Datatype vtypes[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(3, vlens, vdisps, vtypes, &VECTYPE);
    MPI_Type_commit(&VECTYPE);

    return VECTYPE;
}

MPI_Datatype BoidDataType(MPI_Datatype VECTYPE)
{
    Boid boid = Boid(Vec3D(1.0f,2.0f,3.0f), Vec3D(4.0f,5.0f,6.0f), 0);
    MPI_Datatype BOIDTYPE;
    int blens[4] = {1, 1, 1, 1};
    MPI_Aint bdisps[4] = {0};
    MPI_Aint bst_addr, baddr;
    MPI_Get_address(&boid.pos, &bst_addr);
    MPI_Get_address(&boid.vel, &baddr);
    bdisps[1] = baddr - bst_addr;
    MPI_Get_address(&boid.mass, &baddr);
    bdisps[2] = baddr - bst_addr;
    MPI_Get_address(&boid.id, &baddr);
    bdisps[3] = baddr - bst_addr;
    MPI_Datatype btypes[4] = {VECTYPE, VECTYPE, MPI_FLOAT, MPI_INT};
    MPI_Type_create_struct(4, blens, bdisps, btypes, &BOIDTYPE);
    MPI_Type_commit(&BOIDTYPE);

    return BOIDTYPE;
}


int main(int argc, char* argv[])
{
    int argcheck = ArgChecks(argc, argv);
    if (argcheck)
        return 1;

    int threads = std::strtol(argv[1], nullptr, 0);
    int num_boids = std::strtol(argv[2], nullptr, 0);
    std::cout << "Program initialised for " << num_boids << " boids using " << threads << " OpenMP threads." << std::endl;

    omp_set_num_threads(threads); // Set OpenMP to use input no of threads

    MPI_Init(&argc, &argv);
    // Get the number of processes
    int numtasks;
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    int numworkers = numtasks - 1;
    // Get the rank of the process
    int taskid;
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    // Set up tag lists
    int num_timesteps = DURATION / DT;
    std::vector<int> BEGIN_TAGS(num_timesteps);
    std::vector<int> WIND_TAGS(num_timesteps);
    std::vector<int> DONE_TAGS(num_timesteps);
    std::iota(BEGIN_TAGS.begin(), BEGIN_TAGS.end(), 0);
    std::iota(WIND_TAGS.begin(), WIND_TAGS.end(), num_timesteps+1);
    std::iota(DONE_TAGS.begin(), DONE_TAGS.end(), 2*(num_timesteps+1));

    // Set up MPI Datatypes
    MPI_Status status;
    MPI_Datatype VECTYPE = VecDataType();
    MPI_Datatype BOIDTYPE = BoidDataType(VECTYPE);

    // Set up output files
    std::ofstream myfile;
    std::ofstream timingfile;

    // Set up factors to determine how many boids are dealt with on each core
    int split_by = num_boids / numworkers;
    int remaining = num_boids % numworkers;
    if (remaining != 0)
    {
        num_boids += (numworkers - remaining); // Force the number of boids to be a multiple of the number of cores used for easy splitting 
        std::cout << "Adding " << num_boids << " boids for ability to even flock splitting." << std::endl;
        split_by = num_boids / numworkers; // Recalculate split_by
    }

    Vec3D wind_force;
    std::vector<double> dist_matrix;
    std::vector<Boid> flock;

    float TIME = 0;
    int iters = 0;  

    auto tstart = std::chrono::high_resolution_clock::now();

    MPI_Barrier(MPI_COMM_WORLD);
    while (TIME < DURATION)
    {
        if (taskid == DIRECTOR)
        {
            // Begin Director Code

            if (TIME == 0)
            {
                if (numworkers > MAXWORKER || numworkers < MINWORKER)
                {
                    std::cout << "The number of MPI Processors must be between " << MINWORKER+1 << " and " << MAXWORKER+1 << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }

                flock.reserve(num_boids);
                // Generate the flock in a vector on the heap
                auto t1 = std::chrono::high_resolution_clock::now();
                flock = GenFlock(num_boids, POS_ULIM, POS_LLIM, VEL_ULIM, VEL_LLIM);
                auto t2 = std::chrono::high_resolution_clock::now();
                auto time = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1);
                std::cout << "Generated flock with " << flock.size() << " boids in " << time.count()/1e6 << " seconds.\n";

                myfile.open ("pos_dat.csv");
                myfile << "Time,Frame,ID,Mass,X,Y,Z,VelX,VelY,VelZ\n";

                timingfile.open ("timing.csv");
                timingfile << "Frame,Neighbour,Simulation\n";

                // Generate a random wind force
                wind_force = RandWindForce();
            }
            
            // Send full flock and wind force to all workers
            for (int w_id = 1; w_id < numworkers+1; w_id++)
            {
                MPI_Send(&flock.front(), num_boids, BOIDTYPE, w_id, BEGIN_TAGS[iters], MPI_COMM_WORLD);
                MPI_Send(&wind_force, 1, VECTYPE, w_id, WIND_TAGS[iters], MPI_COMM_WORLD);
            }
            
            flock.clear(); // Clear all flock data to ensure new data is entirely from workers
            flock.reserve(num_boids);
            flock.resize(num_boids);

            for (int w_id = 1; w_id < numworkers+1; w_id++)
            {
                int offset_start = (w_id-1) * split_by;
                MPI_Recv(&flock[offset_start], split_by, BOIDTYPE, w_id, DONE_TAGS[iters], MPI_COMM_WORLD, &status);
            }

            for (int num = 0; num < num_boids; num++)
                myfile << TIME << "," << iters << "," << flock[num].mass << "," << flock[num].id << "," << flock[num].pos << "," << flock[num].vel << "\n";

            wind_force = WindEvo(wind_force);
            ProgBar(TIME);

            // End Director code
        }

        else
        {
            // Begin worker code
            
            flock.resize(num_boids);
            Vec3D wind_force;

            int flock_start = (taskid-1) * split_by;
            int flock_end = flock_start + split_by;

            //
            MPI_Recv(&flock.front(), num_boids, BOIDTYPE, DIRECTOR, BEGIN_TAGS[iters], MPI_COMM_WORLD, &status);
            MPI_Recv(&wind_force, 1, VECTYPE, DIRECTOR, WIND_TAGS[iters], MPI_COMM_WORLD, &status);

            dist_matrix.resize(split_by * num_boids);
            dist_matrix = FindDists(flock, flock_start, flock_end, split_by);
            flock = Simulate(flock, flock_start, flock_end, dist_matrix, wind_force, TIME);

            MPI_Send(&flock[flock_start], split_by, BOIDTYPE, DIRECTOR, DONE_TAGS[iters], MPI_COMM_WORLD);

            flock.clear();

            // End Worker Code
        }

        MPI_Barrier(MPI_COMM_WORLD);

        TIME += DT;
        iters += 1;
    }
    

    // Release the memory for custom datatypes
    MPI_Type_free(&BOIDTYPE);
    MPI_Type_free(&VECTYPE);
    // Finalize the MPI environment.
    MPI_Finalize();

    if (taskid == DIRECTOR)
    {
        ProgBar(TIME);
        myfile.close();
        timingfile.close();
        auto tend = std::chrono::high_resolution_clock::now();
        auto timetaken = std::chrono::duration_cast<std::chrono::microseconds>(tend-tstart);
        std::cout << "\nCompleted Simulation for " << DURATION/DT << " timesteps in " << timetaken.count()/1e6 << " seconds.\n";
    }
    

    return 0;
}