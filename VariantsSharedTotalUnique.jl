
# ~~~~~~~~~~~~~~~~~~~~~~~~~ load libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using StatsBase
using StatsPlots
using Statistics
using Measures
# ~~~~~~~~~~~~~~~~~~~~~~ function to read output files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
function read_summed_files(files::Vector{String})
    results_simulation = Dict{Int, DataFrame}()
    
    for (i, file) in enumerate(files)
        if isfile(file)
            df = CSV.read(file, DataFrame)
            results_simulation[i] = df
        else
            println("File not found: $file")
        end
    end
    
    return results_simulation
end
# ~~~~~~~~~~~~~~~~~~~~~~ parameters with which the simulations were run ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
N = 10^3
u = 10^-4
nu = 10^-4
nu_ = round(nu, digits = 6)
m_values = [0.0001, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004] 
nGen = 500000
reps = 5000
begin # theta values 
    θ1 = 0.001
    θ2 = 0.002
    θ3 = 0.003
    θ4 = 0.004
    θ5 = 0.005 
    θ6 = 0.006
    θ7 = 0.007
    θ8 = 0.008
    θ9 = 0.009
    θ10 = 0.01
    θ11 = 0.011 
    θ12 = 0.012  
    θ13 = 0.013  
    θ14 = 0.014 
    θ15 = 0.015 
    θ16 = 0.016  
    θ17 = 0.017  
    θ18 = 0.018  
    θ19 = 0.019  
    θ20 = 0.02 
end
# ~~~~~~~~~~~~~~~~~~~~~~ specify directories ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
dir = "/Documents/MIGRATION/OUTPUT/Variants"
files = ["/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.0001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.0005_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.0015_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.002_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.0025_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.003_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.0035_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv",
         "/Documents/MIGRATION/OUTPUT/Variants/sum_VariantsSummary_N_1000_m_0.004_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02_2.csv"]
# ~~~~~~~~~~~~~~~~~~~~~~ read the files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
results = read_summed_files(files)
for i in eachindex(results)
    results[i] = results[i] ./ reps
end
# ~~~~~~~~~~~~~~~~~~~~~~ functions and colors for plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~ 

function prepare_dicts(shVar_columns, tVar_columns, nVar_columns)
    shared_dict = Dict()
    total_dict = Dict()
    unique_dict = Dict()
    for (i, _) in enumerate(m_values)
        results_i = results[i]
        shared_dict[i] = results_i[:, shVar_columns] |> DataFrame
        total_dict[i] = results_i[:, tVar_columns] |> DataFrame

        unique_df = DataFrame()
        for (nVar_col, shVar_col) in zip(nVar_columns, shVar_columns)
            unique_df[!, Symbol(string(nVar_col))] = (results_i[!, nVar_col] .- results_i[!, shVar_col])
        end
        unique_dict[i] = unique_df
    end
    return unique_dict, shared_dict, total_dict
end
function custom_spectral(n) # Removes yellow color in the middle of the pallette
    full = cgrad(:Spectral, 100, categorical=false)
    left = full.colors[1:40]      
    right = full.colors[61:end]    
    truncated_colors = vcat(left, right)
    return cgrad(truncated_colors, n, categorical=true)
end
function make_plot_set(nVar_cols, shVar_cols, tVar_cols, labels, unique_dict, shared_dict, total_dict)
    plots = []
    for (nVar_col, shVar_col, tVar_col, theta_label) in zip(nVar_cols, shVar_cols, tVar_cols, labels)
        p = plot()
        nVar_values = [unique_dict[m][!, nVar_col] for m in 1:length(m_values)]
        shVar_values = [shared_dict[m][!, shVar_col] for m in 1:length(m_values)]
        tVar_values = [total_dict[m][!, tVar_col] for m in 1:length(m_values)]

        nVar_means = [mean(vals) for vals in nVar_values]
        shVar_means = [mean(vals) for vals in shVar_values]
        tVar_means = [mean(vals) for vals in tVar_values]

        plot!(m_values, nVar_means, linewidth=2.5, color=:grey, label=false, tickfontsize=10, xrotation=45, 
              legend=false, xticks=(m_values, string.(m_values)))
        plot!(m_values, shVar_means, linewidth=2.5, color=:grey, label=false)
        plot!(m_values, tVar_means, linewidth=2.5, color=:grey, label=false)

        scatter!(m_values, nVar_means, marker=:utriangle, markersize=8, markerstrokewidth=0, color=color_map, label=false)
        scatter!(m_values, shVar_means, marker=:circle, markersize=8, markerstrokewidth=0, color=color_map, label=false)
        scatter!(m_values, tVar_means, marker=:diamond, markersize=8, markerstrokewidth=0, color=color_map, label=false)

        push!(plots, p)
    end
    return plots
end
my_colors = custom_spectral(9)  

# the figure will be divided in two rows and three columns 
# specify columns to analyse 
theta_labels_row1 = ["Gen", "Unb"]
theta_labels_row2 = [θ1, θ15, θ20]
shVar_columns_row1 = [:shVarGen, :shVarUnb]
tVar_columns_row1 = [:tVarGen, :tVarUnb]
nVar_columns_row1 = [:nVarP1Gen, :nVarP1Unb]
shVar_columns_row2 = [:shVarCon1, :shVarCon15, :shVarCon20]
tVar_columns_row2 = [:tVarCon1, :tVarCon15, :tVarCon20]
nVar_columns_row2 = [:nVarP1Con1, :nVarP1Con15, :nVarP1Con20]
# extract columns from results and pass to dictionaries 
unique_dict1, shared_dict1, total_dict1 = prepare_dicts(shVar_columns_row1, tVar_columns_row1, nVar_columns_row1)
unique_dict2, shared_dict2, total_dict2 = prepare_dicts(shVar_columns_row2, tVar_columns_row2, nVar_columns_row2)
# make plots 
row1_plots = make_plot_set(nVar_columns_row1, shVar_columns_row1, tVar_columns_row1, theta_labels_row1, unique_dict1, shared_dict1, total_dict1)
row2_plots = make_plot_set(nVar_columns_row2, shVar_columns_row2, tVar_columns_row2, theta_labels_row2, unique_dict2, shared_dict2, total_dict2)
empty_spot = plot(title="", framestyle=:none, grid=false, axis=false)
row1_plots_with_blank = vcat(row1_plots, [empty_spot])  

p = plot(row1_plots_with_blank..., row2_plots...; layout=(2, 3), size=(1800, 1000), bottom_margin=10mm,  left_margin=15mm, right_margin = 10mm)

# save figure
# savefig(p, "$dir/together_total_shared_unique.png")
# savefig(p, "$dir/together_total_shared_unique.svg")
