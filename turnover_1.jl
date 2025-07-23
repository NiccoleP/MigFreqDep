# ~~~~~~~~~~~~~~~~~~~~~~ libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using CategoricalArrays
using StatsBase
using StatsPlots
using Statistics
ENV["GKSwstype"] = "nul"
# ~~~~~~~~~~~~~~~~~~~~~~ functions ~~~~~~~~~~~~~~~~~~~~~~~~~~ 

function calculate_turnover_rate(column_data) # calculates turnover 
    transitions = sum(column_data[i] != column_data[i-1] for i in 2:length(column_data))
    turnover_rate = transitions / length(column_data)
    return turnover_rate
end

function get_turnover(dir_path::AbstractString) # reads files and applies function calculate_turnover_rate
    colnames = [:nVarP1Gen, :nVarP2Gen, :shVarGen, :abVarP1Gen, :fabVarP1Gen, :abVarP2Gen, :fabVarP2Gen, :tVarGen,
                :nVarP1Unb, :nVarP2Unb, :shVarUnb, :abVarP1Unb, :fabVarP1Unb, :abVarP2Unb, :fabVarP2Unb, :tVarUnb,
                :nVarP1Con1, :nVarP2Con1, :shVarCon1, :abVarP1Con1, :fabVarP1Con1, :abVarP2Con1, :fabVarP2Con1, :tVarCon1,
                :nVarP1Con2, :nVarP2Con2, :shVarCon2, :abVarP1Con2, :fabVarP1Con2, :abVarP2Con2, :fabVarP2Con2, :tVarCon2,
                :nVarP1Con3, :nVarP2Con3, :shVarCon3, :abVarP1Con3, :fabVarP1Con3, :abVarP2Con3, :fabVarP2Con3, :tVarCon3,
                :nVarP1Con4, :nVarP2Con4, :shVarCon4, :abVarP1Con4, :fabVarP1Con4, :abVarP2Con4, :fabVarP2Con4, :tVarCon4,
                :nVarP1Con5, :nVarP2Con5, :shVarCon5, :abVarP1Con5, :fabVarP1Con5, :abVarP2Con5, :fabVarP2Con5, :tVarCon5,
                :nVarP1Con6, :nVarP2Con6, :shVarCon6, :abVarP1Con6, :fabVarP1Con6, :abVarP2Con6, :fabVarP2Con6, :tVarCon6,
                :nVarP1Con7, :nVarP2Con7, :shVarCon7, :abVarP1Con7, :fabVarP1Con7, :abVarP2Con7, :fabVarP2Con7, :tVarCon7,
                :nVarP1Con8, :nVarP2Con8, :shVarCon8, :abVarP1Con8, :fabVarP1Con8, :abVarP2Con8, :fabVarP2Con8, :tVarCon8,
                :nVarP1Con9, :nVarP2Con9, :shVarCon9, :abVarP1Con9, :fabVarP1Con9, :abVarP2Con9, :fabVarP2Con9, :tVarCon9,
                :nVarP1Con10, :nVarP2Con10, :shVarCon10, :abVarP1Con10, :fabVarP1Con10, :abVarP2Con10, :fabVarP2Con10, :tVarCon10,
                :nVarP1Con11, :nVarP2Con11, :shVarCon11, :abVarP1Con11, :fabVarP1Con11, :abVarP2Con11, :fabVarP2Con11, :tVarCon11,
                :nVarP1Con12, :nVarP2Con12, :shVarCon12, :abVarP1Con12, :fabVarP1Con12, :abVarP2Con12, :fabVarP2Con12, :tVarCon12,
                :nVarP1Con13, :nVarP2Con13, :shVarCon13, :abVarP1Con13, :fabVarP1Con13, :abVarP2Con13, :fabVarP2Con13, :tVarCon13,
                :nVarP1Con14, :nVarP2Con14, :shVarCon14, :abVarP1Con14, :fabVarP1Con14, :abVarP2Con14, :fabVarP2Con14, :tVarCon14,
                :nVarP1Con15, :nVarP2Con15, :shVarCon15, :abVarP1Con15, :fabVarP1Con15, :abVarP2Con15, :fabVarP2Con15, :tVarCon15,
                :nVarP1Con16, :nVarP2Con16, :shVarCon16, :abVarP1Con16, :fabVarP1Con16, :abVarP2Con16, :fabVarP2Con16, :tVarCon16,
                :nVarP1Con17, :nVarP2Con17, :shVarCon17, :abVarP1Con17, :fabVarP1Con17, :abVarP2Con17, :fabVarP2Con17, :tVarCon17,
                :nVarP1Con18, :nVarP2Con18, :shVarCon18, :abVarP1Con18, :fabVarP1Con18, :abVarP2Con18, :fabVarP2Con18, :tVarCon18,
                :nVarP1Con19, :nVarP2Con19, :shVarCon19, :abVarP1Con19, :fabVarP1Con19, :abVarP2Con19, :fabVarP2Con19, :tVarCon19,
                :nVarP1Con20, :nVarP2Con20, :shVarCon20, :abVarP1Con20, :fabVarP1Con20, :abVarP2Con20, :fabVarP2Con20, :tVarCon20]

    files = readdir(dir_path, join=true)
    varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
    selected_files = varsum_files[1:min(1000, length(varsum_files))] # reads 1000 files 

    data = [rename!(CSV.read(f, DataFrame, header=false), colnames) for f in selected_files]

    all = Vector{Vector{Float64}}()
    for df in data
        row = Float64[]
        for col in names(df)
            if startswith(string(col), "abVarP1") && string(col) != "abVarP1Gen" # get columns with the most abundant variant per trait
                rate = calculate_turnover_rate(df[!, col])
                push!(row, rate)
            end
        end
        push!(all, row)
    end

    return hcat(all...)'  
end

turnover1 = get_turnover("RESULTS/N_1000_m_0.002_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20")
turnover2 = get_turnover("RESULTS/N_1000_m_0.002_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21")

x_values = ["U","θ1","θ2","θ3","θ4","θ5","θ6","θ7","θ8","θ9","θ10","θ11","θ12","θ13","θ14","θ15","θ16","θ17","θ18","θ19","θ20"] # for colnames 

df1 = DataFrame(turnover1, x_values)
CSV.write("df_N_1000_m_0.002_u_0.0001_nu_0.0001_nGen_500000_reps_1000_θsweep_25-05-20.csv", df1)

df2 = DataFrame(turnover2, x_values)
CSV.write("df_N_1000_m_0.002_u_0.0001_nu_0.0005_nGen_500000_reps_1000_θsweep_25-05-21.csv", df2)
