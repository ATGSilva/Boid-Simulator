build: BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp
	@mkdir -p ../../timing
	@mkdir -p ../../boid_data
	@echo Building MPI Boid Simulator...
	@mpicxx $(GCCFLAGS) BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp -o ../../BoidSimMPI
	@sleep 0.1
	@echo "Done."