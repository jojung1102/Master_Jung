# Bending System Analysis

This repository contains the files used for the dynamical and bifurcation analysis. The work involves both 2D and 4D models, numerical simulations, bifurcation diagrams, and trajectory visualizations. Each set of scripts and notebooks is organized by its function in the analysis.

## File Descriptions

### Dynamical System Definitions

- **Bending_system_2D.m**, **Bending_system_4D.m**  
  These MATLAB files define the dynamical equations needed for the bifurcation analysis using MatCont.

### Bifurcation Analysis

- **Bending_system_analysis_2D.m**, **Bending_system_analysis_4D.m**  
  These scripts perform the bifurcation analysis on the 2D and 4D systems, respectively, using the MatCont tool.

### ODE Simulation and Scoring

- **Dynamic_system.m**  
  This file contains the core implementation for running the ODE simulations of the system, including the appropriate time scales.
- **score_main.m**, **score_main_NPC.m**  
  These scripts use the dynamic system module to simulate the ODEs and calculate scores for the results presented in the main text and appendix.

### Trajectory Visualization

- **trajectories_calculation.m**  
  This script calculates example trajectories that are plotted in the main section to demonstrate the dynamic behavior of the system.

### Bifurcation and Flow Visualization

- **bifurcation_flow.ipynb**, **bifurcation_flow_particles.ipynb**  
  These Jupyter notebooks visualize the bifurcation diagrams, the phase flow of the system, and the differences in protein equilibrium for curved proteins compared to uncurved ones.

### Equilibrium Analysis

- **equilibrium_numbers.nb**  
  A Mathematica notebook that calculates explicit expressions for the equilibrium numbers of uncurved proteins.
- **number_equilib.ipynb**  
  A Jupyter notebook that plots the number of equilibrium solutions across different parameters based on the analytical expressions.

### Geometry Visualization

- **Shapes.nb**, **Shape.ipynb**  
  These notebooks visualize the system’s geometry and illustrate the spatial and structural properties related to curvature.

## Notes on Geometry Warnings

Several warnings related to the geometry of the system are present and are discussed in the main text. The reasons for these warnings are marked in the code to provide clarity.
