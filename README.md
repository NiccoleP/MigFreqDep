# Interaction of migration and frequency dependence 
This repository contains the scripts to reproduce all results and plots in 

**citation of paper**

- Simulation of migration between populations: MigrationSimulation.jl  (Julia version 1.10.2 (2024-03-01))

Plotting code for all figures in the manuscript (Julia version 1.11.3 (2025-01-21))

- Figure 1: EquilibriumPhasePlot.jl 
- Figure 2: heatmapJ1vsTheta.jl 
- Figure 3: heatmapTimevsTheta.jl 
- Figure 4: VariantsSharedTotalUnique.jl

Figures in the Supplementary were produced with the same code as the figures in the main text, the difference was the input used (which comes from runs with different parameter setting). 

Plotting code for Supplementary figures: (Julia version 1.10.2 (2024-03-01))

- J0 heatmap -> included in heatmapJ1vsTheta.jl
- Turnover rates: turnover_1.jl reads individual replicate files and outputs a file containing turnover rates (this is done because the script is run in a server and the output files are hgue, therefore it is easier to read and summarize the individual files inside the server and only export the turnover rates); turnover_2.jl reads the output from turnover_1.jl and plots it. 
- Frequency of the most common variant: FreqMostAbVariant.jl

