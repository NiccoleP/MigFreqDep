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
files_ibd = ["/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.0001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.0005_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.001_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.0015_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.002_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.0025_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.003_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv", 
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.0035_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv",
             "/Users/niccole/Documents/MIGRATION/OUTPUT/Variants/sum_IBD_N_1000_m_0.004_u_0.0001_nu_0.0001_reps_5000_thetaθ10.001toθ110.02.csv"]
# ~~~~~~~~~~~~~~~~~~~~~~ read the files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
results_ibd = read_summed_files(files_ibd)
for i in eachindex(results_ibd)
    results_ibd[i] = results_ibd[i] ./ reps
end
#  make dictionary of mean J1 values 
J1_values = Dict()
for i in 1:length(m_values)
    m_ = round(m_values[i], digits = 4)
    rsims = results_ibd[i][49000:50000,:] # get mean over the last 10,000 generations 
    J1_values[m_] = [mean(rsims.J1Gen),
                    mean(rsims.J1Unb),
                    mean(rsims.J1Con1), 
                    mean(rsims.J1Con2),  
                    mean(rsims.J1Con3),
                    mean(rsims.J1Con4),
                    mean(rsims.J1Con5),
                    mean(rsims.J1Con6), 
                    mean(rsims.J1Con7),  
                    mean(rsims.J1Con8),
                    mean(rsims.J1Con9),
                    mean(rsims.J1Con10),
                    mean(rsims.J1Con11),
                    mean(rsims.J1Con12),
                    mean(rsims.J1Con13),
                    mean(rsims.J1Con14),
                    mean(rsims.J1Con15),
                    mean(rsims.J1Con16),
                    mean(rsims.J1Con17),
                    mean(rsims.J1Con18),
                    mean(rsims.J1Con19),
                    mean(rsims.J1Con20)]
end
# sort and transform data 
sorted_m_values = sort(collect(keys(J1_values)))  
J1_means = hcat([J1_values[m] for m in sorted_m_values]...) 
# specify x and y axis 
xlabel = string.(theta_values)
ylabel = string.(m_values)
color_palette = cgrad(:Spectral,10000)
# plot heatmap 
hmp = heatmap(transpose(J1_means),
        xticks=(1:length(theta_values), xlabel),  
        yticks=(1:length(m_values), ylabel),    
        fill_z=transpose(J1_means), 
        aspect_ratio=:equal,
        xlabel = "\$\\theta\$",
        ylabel = "\$m\$",
        colorbar_title = "\$Similarity \\ between \\ populations \\ \\ J_1 \$",
        color=color_palette,
        tickfont=font(6), xrotation=45, dpi= 1000)

# markers for parameters that didn't reach equilibrium
infiniteTime_markers = Dict(0.0001 => [θ4, θ5, θ6, θ7, θ8], 
                        0.0005 => [θ7, θ8, θ9, θ10, θ11],    
                        0.001 => [θ10, θ11, θ12, θ13, θ14],
                        0.0015 => [θ12, θ13, θ14, θ15, θ16, θ17],
                        0.002 => [θ15, θ16, θ17, θ18],
                        0.0025 => [θ17, θ18, θ19, θ20],
                        0.003 => [θ19, θ20])
infiniteEquilibrium_markers = Dict(0.0001 => [θ9, θ10, θ11], 
                                0.0005 => [θ12, θ13],    
                                0.001 => [θ15, θ16],
                                0.0015 => [θ18, θ19],
                                0.002 => [θ19, θ20])

# add markers to plot 
for (m_val, θ_list) in infiniteTime_markers
    m_index = findfirst(x -> x == m_val, sorted_m_values)
    if m_index !== nothing
        θ_indices_mapped = findall(x -> x in θ_list, theta_values)
        for θ_idx in θ_indices_mapped
            scatter!([θ_idx], [m_index], 
            markershape=:hexagon, 
            color=RGB(40/255, 110/255, 156/255), 
            markersize=3, 
            markerstrokewidth=0, legend = false) 
            end
        end
end              
for (m_val, θ_list) in infiniteEquilibrium_markers
    m_index = findfirst(x -> x == m_val, sorted_m_values)
        if m_index !== nothing
            θ_indices_mapped = findall(x -> x in θ_list, theta_values)
            for θ_idx in θ_indices_mapped
                scatter!([θ_idx], [m_index], 
                markershape=:diamond, 
                color=RGB(255/255, 188/255, 66/255), 
                markersize=3, 
                markerstrokewidth=0, legend = false) 
            end
        end
end 

# save plot 
savefig(hmp, "$dir/hmpJ1_m_vs_theta.png")
savefig(hmp, "$dir/hmpJ1_m_vs_theta.svg")

# Similarly, for J0 heatmap:

# make dictionary for J0 values 
J0_values = Dict()
for i in 1:length(m_values)
    m_ = round(m_values[i], digits = 4)
    rsims = results_ibd[i][49000:50000,:] # get mean over the last 10,000 generations 
    J0_values[m_] = [mean(rsims.J0Gen1),
                    mean(rsims.J0Unb1),
                    mean(rsims.J0Con1P1), 
                    mean(rsims.J0Con2P1),  
                    mean(rsims.J0Con3P1),
                    mean(rsims.J0Con4P1),
                    mean(rsims.J0Con5P1),
                    mean(rsims.J0Con6P1), 
                    mean(rsims.J0Con7P1),  
                    mean(rsims.J0Con8P1),
                    mean(rsims.J0Con9P1),
                    mean(rsims.J0Con10P1),
                    mean(rsims.J0Con11P1),
                    mean(rsims.J0Con12P1),
                    mean(rsims.J0Con13P1),
                    mean(rsims.J0Con14P1),
                    mean(rsims.J0Con15P1),
                    mean(rsims.J0Con16P1),
                    mean(rsims.J0Con17P1),
                    mean(rsims.J0Con18P1),
                    mean(rsims.J0Con19P1),
                    mean(rsims.J0Con20P1)]
end

theta_values = ["G","U",θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20]
sorted_m_values = sort(collect(keys(J0_values)))  
J0_means = hcat([J0_values[m] for m in sorted_m_values]...) 

xlabel = string.(theta_values)
ylabel = string.(m_values)

color_palette = cgrad(:Spectral,10000)
hmp = heatmap(transpose(J0_means),
        xticks=(1:length(theta_values), xlabel),  
        yticks=(1:length(m_values), ylabel),    
        fill_z=transpose(J0_means), 
        aspect_ratio=:equal,
        xlabel = "\$\\theta\$",
        ylabel = "\$m\$",
        colorbar_title = "\$Similarity \\ within \\ populations \\ \\ J_0 \$",
        color=color_palette,
        tickfont=font(6), xrotation=45, dpi= 1000)

savefig(hmp, "$dir/hmpJ0_m_vs_theta.png")
savefig(hmp, "$dir/hmpJ0_m_vs_theta.svg")



# Calculating I 

all_I = Vector{Vector{Float64}}(undef, length(m_values))
all_D = Vector{Vector{Float64}}(undef, length(m_values))

for i in 1:length(m_values)
        m_ = m_values[i]
        rsims = results_ibd[i][49000:50000,:]
        J0_values= [mean(rsims.J0Gen1), # get J0 values 
                    mean(rsims.J0Unb1),
                    mean(rsims.J0Con1P1), 
                    mean(rsims.J0Con2P1),  
                    mean(rsims.J0Con3P1),
                    mean(rsims.J0Con4P1),
                    mean(rsims.J0Con5P1),
                    mean(rsims.J0Con6P1), 
                    mean(rsims.J0Con7P1),  
                    mean(rsims.J0Con8P1),
                    mean(rsims.J0Con9P1),
                    mean(rsims.J0Con10P1),
                    mean(rsims.J0Con11P1),
                    mean(rsims.J0Con12P1),
                    mean(rsims.J0Con13P1),
                    mean(rsims.J0Con14P1),
                    mean(rsims.J0Con15P1),
                    mean(rsims.J0Con16P1),
                    mean(rsims.J0Con17P1),
                    mean(rsims.J0Con18P1),
                    mean(rsims.J0Con19P1),
                    mean(rsims.J0Con20P1)]

        J1_values= [mean(rsims.J1Gen), # get J1 values
                    mean(rsims.J1Unb),
                    mean(rsims.J1Con1), 
                    mean(rsims.J1Con2),  
                    mean(rsims.J1Con3),
                    mean(rsims.J1Con4),
                    mean(rsims.J1Con5),
                    mean(rsims.J1Con6), 
                    mean(rsims.J1Con7),  
                    mean(rsims.J1Con8),
                    mean(rsims.J1Con9),
                    mean(rsims.J1Con10),
                    mean(rsims.J1Con11),
                    mean(rsims.J1Con12),
                    mean(rsims.J1Con13),
                    mean(rsims.J1Con14),
                    mean(rsims.J1Con15),
                    mean(rsims.J1Con16),
                    mean(rsims.J1Con17),
                    mean(rsims.J1Con18),
                    mean(rsims.J1Con19),
                    mean(rsims.J1Con20)]
        
        I = J1_values ./ J0_values # calculate ratio 
        #I = J1_values ./ (J0_values).^1/2
        D = -log.(I) # Nei's distance
        all_I[i] = I
        all_D[i] = D
end

D_matrix = hcat(all_D...)'  # rows = m_values, columns = traits (transpose to match)
xlabel = string.(theta_values)
ylabel = string.(m_values)
color_palette = cgrad(:Spectral,10000)
d_distance = heatmap(1:length(theta_values),        
        1:length(m_values),            
        D_matrix,                     
        xticks=(1:length(theta_values), xlabel),
        yticks=(1:length(m_values), ylabel),
        xlabel = "\$traits\$",
        ylabel = "\$m\$",
        colorbar_title = "D",
        color = color_palette,
        aspect_ratio = :equal,  xrotation=45, dpi = 2000)
savefig(d_distance,"$dir/d_distance.svg")
savefig(d_distance,"$dir/d_distance.png")


# D_gen = [d[1] for d in all_D]   # index 1 = Gen
# D_unb = [d[2] for d in all_D]   # index 2 = Unb
# D_con1 = [d[3] for d in all_D]
# D_con2 = [d[4] for d in all_D]
# D_con3 = [d[5] for d in all_D]
# D_con4 = [d[6] for d in all_D]
# D_con5 = [d[7] for d in all_D]
# D_con6 = [d[8] for d in all_D]
# D_con7 = [d[9] for d in all_D]
# D_con8 = [d[10] for d in all_D]
# D_con9 = [d[11] for d in all_D]
# D_con10 = [d[12] for d in all_D]
# D_con11 = [d[13] for d in all_D]
