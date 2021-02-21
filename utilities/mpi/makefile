animation: plot_from_data.py
	@mkdir -p ../../animations
	@echo "Output will be in animations directory."
	@python3 plot_from_data.py
	@echo "Done."

dashboard: plot_dashboard.py FlockingCoeff.py
	@echo "WARNING: This will take a very long time for large N systems (N > 3000)."
	@mkdir -p ../../animations
	@echo "Output will be in animations directory."
	@echo "Calculating flocking coefficients..."
	@python FlockingCoeff.py
	@echo "Producing dashboard animation..."
	@python plot_dashboard.py
	@echo "Done."