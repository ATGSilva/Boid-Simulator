build: BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp Vec3D.cpp
	@mkdir -p ../../timing
	@mkdir -p ../../boid_data
	@echo Building OpenMP Boid Simulator...
	@g++ $(GCCFLAGS) BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp Vec3D.cpp -o ../../BoidSimOMP
	@sleep 0.1
	@echo "Done."