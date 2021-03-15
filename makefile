.DEFAULT_GOAL := help
GCCFLAGS = -fopenmp -march=native -mtune=native -O3 -std=c++11 -lm -g
NP = 5
PPN = 2
OMP_THR = 10
NUM_BOIDS = 1000
export

all_omp: build_omp run_omp animation_omp dashboard_omp

all_mpi: build_mpi run_mpi animation_mpi dashboard_mpi

build_omp:
	@+$(MAKE) build -C src/omp -s

build_mpi:
	@+$(MAKE) build -C src/mpi -s

run_omp: BoidSimOMP
	@./BoidSimOMP $(OMP_THR) $(NUM_BOIDS)

run_mpi: BoidSimMPI
	@mpirun -np $(NP) ./BoidSimMPI $(PPN) $(NUM_BOIDS)

clean:
	@echo "Cleaning binary executable, results and animation files."
	@rm -f BoidSimMPI
	@rm -f BoidSimOMP
	@rm -rf boid_data
	@rm -rf timing
	@rm -rf animations
	@echo "Done."

animation_omp: 
	@$(MAKE) animation -C utilities/omp -s

animation_mpi: 
	@$(MAKE) animation -C utilities/mpi -s

dashboard_omp: 
	@+$(MAKE) dashboard -C utilities/omp -s

dashboard_mpi: 
	@+$(MAKE) dashboard -C utilities/mpi -s
	
help:
	@echo "Simulation settings can be changed by editing the 'Settings.h' file, see Readme.txt for more info."
	@echo "\nTo build an OpenMP program use 'make build_omp'."
	@echo "To build an MPI program use 'make build_mpi'."
	@echo "\nTo override the number of boids to simulate use 'make run_<omp/mpi> NUM_BOIDS=<Desired Number>'."
	@echo "To override the number of threads in an OpenMP run use 'make run_omp OMP_THR=<Desired Number>'."
	@echo "To override the number of nodes and threads per node in an MPI run use 'make run_mpi NP=<Desired Nodes> PPN=<Desired Threads per Node>'."
	@echo "\nWARNING: using 'make clean' will remove the current build executable as well as all generated data, animations and dashboards."
	@echo "\nFor help with animation and dashboard creation see Readme.txt."