.DEFAULT_GOAL := build

GCCFLAGS = -fopenmp -march=native -O3

build: BoidSimMPIOMP1.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp
	mpicxx $(GCCFLAGS) BoidSimMPIOMP1.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim

run:
	@echo "To run this file use command 'mpirun -np <Num_MPI_procs> ./BoidSim <Num_OMP_threads> <Num_Boids>'"

clean:
	rm BoidSim