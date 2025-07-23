# ~~~~~~~~~~~~~~~~~~~~~~ libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using StatsBase
using StatsPlots
using Statistics
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
begin
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
    theta_values = ["G","U",θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20]
end
# ~~~~~~~~~~~~~~~~~~~~~~ specify directories ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
dir = "/Documents/MIGRATION/OUTPUT/" # change path to your local directory 
files = ["/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.0001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.0005_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.0015_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.002_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.0025_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.003_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.0035_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
         "/Documents/MIGRATION/OUTPUT/sum_VariantsSummary_N_1000_m_0.004_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv"]
# ~~~~~~~~~~~~~~~~~~~~~~ read the files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
results = read_summed_files(files)
for i in eachindex(results)
    results[i] = results[i] ./ reps
end

#~~~~~~~~~~~~~~~~~~~~~~     Calculating time to equilibrium functions      ~~~~~~~~~~~~~~~~~~~~~~

# Calculates the average value of data from start_gen to the end and finds the first index where the data exceeds this average.
function calculate_limit(data,start_gen)
         limitValue = mean(data[start_gen:end])
         first_exceed_index = findfirst(x -> x > limitValue, data)
         return (limitValue, first_exceed_index)
end

# This function creates a dictionary to store the equilibrium values for all J1 columns in a dataframe. 
# It applies function 'calculate_limit' to each J1 column. It returns the dictionary with the equilibrium values and a plot of the data.
# If you want to have the plot output, simply uncomment the code inside the function 
# It returns the equilibrium generation and a plot
function get_equilibrium_times_j1(results::DataFrame,start_gen)
    # plots = []
    equilibrium = Dict{String, Union{Int, Float64}}()  
    for col in names(results)  
        if occursin("J1", col) 
            data = results[:, col] 
            limitValue, generation = calculate_limit(data,start_gen)
            equilibrium[col] = generation
            # p = plot(data, title=col, label = false, xtickfontsize=6)
            # if !isnothing(equilibrium)
            #     vline!(p, [generation], label = false)
            #     hline!(p, [limitValue], label = false)
            # end
            # push!(plots, p)  
        end
    end
    return equilibrium #, plots 
end

# Sets equilibrium values for specified keys to NaN in the equilibrium dictionary
function markNothing!(keys_to_modify, equilibrium_times)
    for key in keys_to_modify
        if haskey(equilibrium_times, key)
            equilibrium_times[key] = NaN
        end
    end
end

equilibrium_times_j1_m1 = get_equilibrium_times_j1(results[1],450000) # if one wants to see the plot add plots_m1 before the equal sign 
equilibrium_times_j1_m2 = get_equilibrium_times_j1(results[2],450000)
equilibrium_times_j1_m3 = get_equilibrium_times_j1(results[3],450000)
equilibrium_times_j1_m4 = get_equilibrium_times_j1(results[4],450000)
equilibrium_times_j1_m5 = get_equilibrium_times_j1(results[5],450000)
equilibrium_times_j1_m6 = get_equilibrium_times_j1(results[6],450000)
equilibrium_times_j1_m7 = get_equilibrium_times_j1(results[7],450000)
equilibrium_times_j1_m8 = get_equilibrium_times_j1(results[8],450000)
equilibrium_times_j1_m9 = get_equilibrium_times_j1(results[9],450000)

# order in which I want the columns  
key_order = ["J1Gen","J1Unb","J1Con1", "J1Con2", "J1Con3", "J1Con4", "J1Con5", "J1Con6", "J1Con7", "J1Con8", "J1Con9", "J1Con10", 
             "J1Con11",  "J1Con12", "J1Con13", "J1Con14", "J1Con15", "J1Con16", "J1Con17", "J1Con18", "J1Con19", "J1Con20"]
# Construct a matrix where each column corresponds to the equilibrium times from a different migration condition (m1 to m9) ordered by key_order.
# horizontally concatenate each set of equilibrium times as new columns
matrix = [get(equilibrium_times_j1_m1, key, NaN) for key in key_order]
matrix = hcat(matrix, [get(equilibrium_times_j1_m2, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m3, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m4, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m5, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m6, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m7, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m8, key, NaN) for key in key_order])
matrix = hcat(matrix, [get(equilibrium_times_j1_m9, key, NaN) for key in key_order])

matrix = matrix'
color_palette = cgrad(:Spectral)
theta_values = ["G","U",θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20]

xlabel = string.(theta_values)
ylabel = string.(m_values)
hmp_time = heatmap(matrix,
        xticks=(1:length(theta_values), xlabel),  
        yticks=(1:length(m_values), ylabel),    
        fill_z=transpose(matrix), 
        aspect_ratio=:equal,
        xlabel = "\$\\theta\$",
        ylabel = "\$m\$",
        colorbar_title = "\$Time \\ to \\ equilibrium \$",
        color=color_palette,
        colorbar_scale=:log10,
        tickfont=font(6), xrotation=45, dpi = 1000)


        for (m_val, θ_list) in infiniteTime_markers # add marker for parameter combinations that never reached equilibrium 
            m_index = findfirst(x -> x == m_val, m_values)
            if m_index !== nothing
                θ_indices_mapped = findall(x -> x in θ_list, theta_values)
                for θ_idx in θ_indices_mapped
                    scatter!([θ_idx], [m_index], 
                        markershape=:hexagon, 
                        color=RGB(255/255, 184/255, 77/255), 
                        markersize=3, 
                        markerstrokewidth=0, legend = false) 
                end
            end
        end        

        for (m_val, θ_list) in infiniteEquilibrium_markers  # add marker for parameter combinations that can be considered at equilibrium 
            m_index = findfirst(x -> x == m_val, m_values)
            if m_index !== nothing
                θ_indices_mapped = findall(x -> x in θ_list, theta_values)
                for θ_idx in θ_indices_mapped
                    scatter!([θ_idx], [m_index], 
                        markershape=:diamond, 
                        color=RGB(176/255, 226/255, 152/255), 
                        markersize=3, 
                        markerstrokewidth=0, legend = false) 
                end
            end
        end

savefig(hmp_time, "$dir/hmp_time_m_vs_theta.png")
savefig(hmp_time, "$dir/hmp_time_m_vs_theta.svg")