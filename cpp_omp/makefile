.DEFAULT_GOAL := build
GCCFLAGS = -fopenmp -march=native -O3
NUM_THREADS = 10
NUM_BOIDS = 1000

all: build run animation

build: BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp
	@mkdir -p timing
	@mkdir -p boid_data
	@echo Building...
	@g++ $(GCCFLAGS) BoidSimOMP.cpp Neighbours.cpp Forces.cpp Boid.cpp Vec3D.cpp -o BoidSim
	@echo "Done."

run: BoidSim
	@./BoidSim $(NUM_THREADS) $(NUM_BOIDS) 

clean:
	@echo "Cleaning binary executable, results and animation files."
	@rm BoidSim
	@rm boid_data/*.csv
	@rm timing/*.csv
	@echo "Done."

animation: plot_from_data.py
	@mkdir -p animations
	@echo "Output will be in animations directory."
	@python3 plot_from_data.py $(NUM_BOIDS) $(NUM_THREADS)
	@echo "Done."