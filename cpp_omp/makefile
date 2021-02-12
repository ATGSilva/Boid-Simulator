GCCFLAGS = -fopenmp -march=native -O3

make: BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp
		g++ $(GCCFLAGS) BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim