# Interplay of migration and frequency dependence 
This repositody contains the scripts to reproduce all results and plots in 

**APA citation of paper**

link to paper 

- Simulation of migration between populations: MigrationSimulation.jl

Plotting code for all figure in the manuscript

- Figure 1: EquilibriumPhasePlot.jl
- Figure 2: heatmapJ1vsTheta.jl
- Figure 3: heatmapTimevsTheta.jl
- Figure 4: VariantsSharedTotalUnique.jl

Figures in the Supplementary were produced with the same code as the figures in the main text, the difference was the input used (which comes from runs with different parameter setting). 

Plotting code for additional Supplementary figures: 

- J0 heatmap -> included in heatmapJ1vsTheta.jl
- Nei's distance heatmap -> included in heatmapJ1vsTheta.jl 
- turnover rates: turnover_1.jl reads individual replicate files and outputs a file containing turnover rates (this is done because the script is run in a server and the output files are hgue, therefore it is easier to read and summarize the individual files inside the server and only export the turnover rates); turnover_2.jl reads the output from turnover_1.jl and plots it. 
- frequency of the most common variant: FreqMostAbVariant.jl

