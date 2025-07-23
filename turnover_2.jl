# ~~~~~~~~~~~~~~~~~~~~~~ libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using StatsBase
using StatsPlots
using Statistics
using CategoricalArrays
# ~~~~~~~~~~~~~~~~~~~~~~ parameters with which the simulations were run ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
N = 10^3
u = 10^-4
nu = 10^-4
nu_ = round(nu, digits = 6)
m_values = [0.0001, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004] 
nGen = 500000
reps = 5000
# ~~~~~~~~~~~~~~~~~~~~~~ specify directories ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
dir = "/Documents/MIGRATION/OUTPUT/turnover/"
files = ["/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.0001_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20.csv",
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.0001_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21.csv",
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.0005_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20.csv",
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.0005_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21.csv",
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.001_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20.csv",
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.001_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21.csv", 
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.002_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20.csv", 
         "/Documents/MIGRATION/OUTPUT/turnover/df_N_1000_m_0.002_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21.csv"]

m_values = ["0.0001", "0.0001", "0.0005", "0.0005", "0.001", "0.001", "0.002", "0.002"]
nu_values = ["0.0001", "0.0005", "0.0001", "0.0005", "0.0001", "0.0005", "0.0001", "0.0005"]
# ~~~~~~~~~~~~~~~~~~~~~~ function to read output files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
m_0001_u_0001 = CSV.read(files[1], DataFrame, header = true) 
m_0001_u_0005 = CSV.read(files[2], DataFrame, header = true)
m_0005_u_0001 = CSV.read(files[3], DataFrame, header = true)
m_0005_u_0005 = CSV.read(files[4], DataFrame, header = true)
m_001_u_0001 = CSV.read(files[5], DataFrame, header = true)
m_001_u_0005 = CSV.read(files[6], DataFrame, header = true)
m_002_u_0001 = CSV.read(files[7], DataFrame, header = true)
m_002_u_0005 = CSV.read(files[8], DataFrame, header = true)
# ~~~~~~~~~~~~~~~~~~~~~~ plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~ 

x_values = ["U","θ1","θ2","θ3","θ4","θ5","θ6","θ7","θ8","θ9","θ10","θ11","θ12","θ13","θ14","θ15","θ16","θ17","θ18","θ19","θ20"]
dfs = [m_0001_u_0001, m_0001_u_0005, m_0005_u_0001, m_0005_u_0005,
       m_001_u_0001, m_001_u_0005, m_002_u_0001, m_002_u_0005]
means_per_df = [mean.(eachcol(df)) for df in dfs]
stds_per_df = [std.(eachcol(df)) for df in dfs]
my_colors = [:blue, :blue, :green, :green, :orange, :orange, :purple, :purple]
my_colors2= palette(:Spectral, 9)
p = plot(dpi = 1000)
for (i, means) in enumerate(means_per_df)
    linestyle = (iseven(i) ? :dash : :solid)  # even = dotted; odd = solid
    label_str = "m = $(m_values[i]), ν = $(nu_values[i])"
    stds = stds_per_df[i]
    plot!(
        p, x_values, means,
        #ribbon = stds, 
        #fillalpha = 0.2,
        xlabel = "\$θ\$",
        ylabel = "\$Turnover\$",
        label = label_str,
        foreground_color_legend = nothing,
        background_color_legend = nothing, 
        lw = 1.5,
        linestyle = linestyle,
        color = my_colors2[i],
        xticks = (1:length(x_values), x_values),
        xrotation = 45)
end

savefig(p,"$dir/turnover_means.png")