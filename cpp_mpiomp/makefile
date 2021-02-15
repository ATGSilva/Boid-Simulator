.DEFAULT_GOAL := build
GCCFLAGS = -fopenmp -march=native -O3
NUM_PROCS = 5
NUM_THREADS = 2
NUM_BOIDS = 1000

all: build run_desktop animation

build: BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp
	@mkdir -p timing
	@mkdir -p boid_data 
	@echo "Building..."
	@mpicxx $(GCCFLAGS) BoidSimMPIOMP.cpp Simulate.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim
	@echo "Done."

run_desktop: BoidSim
	@mpirun -np $(NUM_PROCS) ./BoidSim $(NUM_THREADS) $(NUM_BOIDS)

clean:
	@echo "Cleaning binary executable, results and animation files."
	@rm BoidSim
	@rm boid_data/*.csv
	@rm timing/*.csv
	@rm animations/*.mp4
	@echo "Done."

animation: plot_from_data.py
	@mkdir -p animations
	@echo "Output will be in animations directory."
	@python3 plot_from_data.py $(NUM_BOIDS) $(NUM_THREADS) $(NUM_PROCS)
	@echo "Done."