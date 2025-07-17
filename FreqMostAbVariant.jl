# ~~~~~~~~~~~~~~~~~~~~~~ libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using StatsBase
using StatsPlots
using Statistics
ENV["GKSwstype"] = "nul"
# ~~~~~~~~~~~~~~~~~~~~~~ specify columns and values for reading and plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~ 

colnames = [:nVarP1Gen, :nVarP2Gen, :shVarGen, :abVarP1Gen, :fabVarP1Gen, :abVarP2Gen, :fabVarP2Gen, :tVarsGen,
            :nVarP1Unb, :nVarP2Unb, :shVarUnb, :abVarP1Unb, :fabVarP1Unb, :abVarP2Unb, :fabVarP2Unb, :tVarsUnb,
            :nVarP1Con1, :nVarP2Con1, :shVarCon1, :abVarP1Con1, :fabVarP1Con1, :abVarP2Con1, :fabVarP2Con1, :tVarsCon1,
            :nVarP1Con2, :nVarP2Con2, :shVarCon2, :abVarP1Con2, :fabVarP1Con2, :abVarP2Con2, :fabVarP2Con2, :tVarsCon2,
            :nVarP1Con3, :nVarP2Con3, :shVarCon3, :abVarP1Con3, :fabVarP1Con3, :abVarP2Con3, :fabVarP2Con3, :tVarsCon3,
            :nVarP1Con4, :nVarP2Con4, :shVarCon4, :abVarP1Con4, :fabVarP1Con4, :abVarP2Con4, :fabVarP2Con4, :tVarsCon4,
            :nVarP1Con5, :nVarP2Con5, :shVarCon5, :abVarP1Con5, :fabVarP1Con5, :abVarP2Con5, :fabVarP2Con5, :tVarsCon5,
            :nVarP1Con6, :nVarP2Con6, :shVarCon6, :abVarP1Con6, :fabVarP1Con6, :abVarP2Con6, :fabVarP2Con6, :tVarsCon6,
            :nVarP1Con7, :nVarP2Con7, :shVarCon7, :abVarP1Con7, :fabVarP1Con7, :abVarP2Con7, :fabVarP2Con7, :tVarsCon7,
            :nVarP1Con8, :nVarP2Con8, :shVarCon8, :abVarP1Con8, :fabVarP1Con8, :abVarP2Con8, :fabVarP2Con8, :tVarsCon8,
            :nVarP1Con9, :nVarP2Con9, :shVarCon9, :abVarP1Con9, :fabVarP1Con9, :abVarP2Con9, :fabVarP2Con9, :tVarsCon9,
            :nVarP1Con10, :nVarP2Con10, :shVarCon10, :abVarP1Con10, :fabVarP1Con10, :abVarP2Con10, :fabVarP2Con10, :tVarsCon10,
            :nVarP1Con11, :nVarP2Con11, :shVarCon11, :abVarP1Con11, :fabVarP1Con11, :abVarP2Con11, :fabVarP2Con11, :tVarsCon11,
            :nVarP1Con12, :nVarP2Con12, :shVarCon12, :abVarP1Con12, :fabVarP1Con12, :abVarP2Con12, :fabVarP2Con12, :tVarsCon12,
            :nVarP1Con13, :nVarP2Con13, :shVarCon13, :abVarP1Con13, :fabVarP1Con13, :abVarP2Con13, :fabVarP2Con13, :tVarsCon13,
            :nVarP1Con14, :nVarP2Con14, :shVarCon14, :abVarP1Con14, :fabVarP1Con14, :abVarP2Con14, :fabVarP2Con14, :tVarsCon14,
            :nVarP1Con15, :nVarP2Con15, :shVarCon15, :abVarP1Con15, :fabVarP1Con15, :abVarP2Con15, :fabVarP2Con15, :tVarsCon15,
            :nVarP1Con16, :nVarP2Con16, :shVarCon16, :abVarP1Con16, :fabVarP1Con16, :abVarP2Con16, :fabVarP2Con16, :tVarsCon16,
            :nVarP1Con17, :nVarP2Con17, :shVarCon17, :abVarP1Con17, :fabVarP1Con17, :abVarP2Con17, :fabVarP2Con17, :tVarsCon17,
            :nVarP1Con18, :nVarP2Con18, :shVarCon18, :abVarP1Con18, :fabVarP1Con18, :abVarP2Con18, :fabVarP2Con18, :tVarsCon18,
            :nVarP1Con19, :nVarP2Con19, :shVarCon19, :abVarP1Con19, :fabVarP1Con19, :abVarP2Con19, :fabVarP2Con19, :tVarsCon19,
            :nVarP1Con20, :nVarP2Con20, :shVarCon20, :abVarP1Con20, :fabVarP1Con20, :abVarP2Con20, :fabVarP2Con20, :tVarsCon20]
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
end
fabVars = [:fabVarP1Gen, :fabVarP1Unb, :fabVarP1Con1, :fabVarP1Con2, :fabVarP1Con3, :fabVarP1Con4, :fabVarP1Con5,
          :fabVarP1Con6, :fabVarP1Con7, :fabVarP1Con8, :fabVarP1Con9, :fabVarP1Con10, :fabVarP1Con11,
          :fabVarP1Con12, :fabVarP1Con13, :fabVarP1Con14, :fabVarP1Con15, :fabVarP1Con16, :fabVarP1Con17,
          :fabVarP1Con18, :fabVarP1Con19, :fabVarP1Con20]
x_values = ["G","U", θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20]

# ~~~~~~~~~~~~~~~~~~~~~~ plot per m rate ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
m = 0.004
files = readdir("RESULTS/N_1000_m_0.004_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-03", join=true) # read directory
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files) # filter files 
selected_files = varsum_files[1:min(100, length(varsum_files))] # specify 100 files 
data = [CSV.read(f, DataFrame, header = false) for f in selected_files] # read CSVs
data = [rename!(df, colnames) for df in data] # rename columns 
fabVar_df = DataFrame()
for df in data                          # make dataframe with all the columns involving the most abundant variants per trait 
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], # plot 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg") # save figure 



m = 0.0035
files = readdir("RESULTS/N_1000_m_0.0035_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-03_3", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")



m = 0.003
files = readdir("RESULTS/N_1000_m_0.003_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-04_3", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")



m = 0.0025
files = readdir("RESULTS/N_1000_m_0.0025_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-04", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")


m = 0.002
files = readdir("RESULTS/N_1000_m_0.002_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-04", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")


m = 0.0015
files = readdir("RESULTS/N_1000_m_0.0015_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-07_3", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")


m = 0.001
files = readdir("RESULTS/N_1000_m_0.001_u_0.0001_nu_0.0001_nGen_500000_reps_100_θsweep_25-04-03", join=true)
varsum_files = filter(f -> occursin("VarSumAll_replicate", f), files)
selected_files = varsum_files[1:min(100, length(varsum_files))]
data = [CSV.read(f, DataFrame, header = false) for f in selected_files]
data = [rename!(df, colnames) for df in data]
fabVar_df = DataFrame()
for df in data
    push!(fabVar_df, df[end, fabVars])
end
means = mean.(eachcol(fabVar_df))

p = violin([fabVar_df[!,col] for col in fabVars], 
        legend=false, 
        ylabel = "\$frequency \\ of \\ most \\ popular \\ variant\$",
        xlabel = "\$\\theta\$",
        line = (0.35, :blue), 
        fill = (0.2, :blue))
x_ticks = 1:length(fabVars), x_values
plot!(p, xticks=x_ticks, xrotation=45)
scatter!(p, means, color=:darkblue, marker=:circle, markerstrokewidth=0, markersize=1.5)
savefig(p, "frequency_most_abundant_variant_m_$(m)_reps_100.svg")