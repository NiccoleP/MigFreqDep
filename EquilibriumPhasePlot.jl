# ~~~~~~~~~~~~~~~~~~~~~~~~~ load libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
using Plots
using DataFrames
using CSV
using StatsBase
using StatsPlots
using Statistics
# ~~~~~~~~~~~~~~~~~~~~~~ function to read output files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
function read_files(files::Vector{String}) # function to read files 
    results_simulation = Dict{Int, DataFrame}()
    for (i, file) in enumerate(files)
        if isfile(file)
            df = CSV.read(file, DataFrame)
            results_simulation[i] = df
        else
            println("File not found: $file")  # error in case it can't find the path you input 
        end
    end
    return results_simulation
end
# ~~~~~~~~~~~~~~~~~~~~~~ parameters with which the simulations were run ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
N = 10^3                                    # population size
u = 10^-4                                   # mutation rate
nu = 10^-4                                  # innovation rate
nu_ = round(nu, digits = 6)                 # modify for figure output
m_values = [0.0001, 0.0005, 0.001, 0.002]   # migration rates
nGen = 500000                               # number of generations simulations ran for
reps = 5000                                 # number of replicates
θ0 = 0                                      # unbiased trait
θ1  = -0.020                                # θ1 to θ8 are anti conformity values 
θ2  = -0.015  
θ3  = -0.010  
θ4  = -0.0075  
θ5  = -0.005  
θ6  = -0.0025  
θ7  = -0.001  
θ8  = -0.0005  
θ9 = 0.0005                                # θ9 to θ20 are conformity values 
θ10 = 0.001
θ11 = 0.002
θ12 = 0.004
θ13 = 0.006
θ14 = 0.008
θ15 = 0.010
θ16 = 0.012
θ17 = 0.014
θ18 = 0.016
θ19 = 0.018
θ20 = 0.020
# ~~~~~~~~~~~~~~~~~~~~~~ specify directories ~~~~~~~~~~~~~~~~~~~~~~~~~~ 
dir = "/MigFreqDep/OUTPUT/EqPlot/" # Change directory to whererever you have your output files 
files = ["/MigFreqDep/OUTPUT/EqPlot/sum_IBD_N_1000_m_0.0001_u_0.0001_nu_0.0001_reps_5000_thetaθ1-0.02toθ110.02.csv",
         "/MigFreqDep/MIGRATION/OUTPUT/EqPlot/sum_IBD_N_1000_m_0.0005_u_0.0001_nu_0.0001_reps_5000_thetaθ1-0.02toθ110.02.csv",
         "/MigFreqDep/MIGRATION/OUTPUT/EqPlot/sum_IBD_N_1000_m_0.001_u_0.0001_nu_0.0001_reps_5000_thetaθ1-0.02toθ110.02.csv",
         "/MigFreqDep/MIGRATION/OUTPUT/EqPlot/sum_IBD_N_1000_m_0.002_u_0.0001_nu_0.0001_reps_5000_thetaθ1-0.02toθ110.02.csv",
         "/MigFreqDep/MIGRATION/OUTPUT/EqPlot/sum_IBD_N_1000_m_0_u_0.0001_nu_0.0001_reps_5000_thetaθ1-0.02toθ110.02.csv"] # last file are runs without migration for the rug part on the bottom of the plot 
# ~~~~~~~~~~~~~~~~~~~~~~ read the files ~~~~~~~~~~~~~~~~~~~~~~~~~~ 

results = read_files(files)  # read files 
for i in eachindex(results)  # divide by number of replicates since the output files have the sum of all the individual run's J0 and J1 values
    results[i] = results[i] ./ reps  
end



# GENES EXPECTATION
results_theoretical = Dict()
exp_values_g = Dict() 
speed_g = Dict()
N = 10^3
u = 10^-4
nu = 10^-4
nu_ = round(nu, digits = 6)
begin # Expectations Nei Feldman 1972: J0_exp, J1_exp, eigenvector, and λ
    for m in m_values
        u =  10^-4 # modify here if you want to get the expectation for other μ values, e.g. 5 * 10^-5
        a = (1 - m)^2 + m^2
        G = 1 - (a * (1 - u)^2) * (2 - 1 / (2 * N)) - (1 - 2 * a) * ((1 - u)^4) * (1 - 1 / (2 * N))
        c = (1 - a + u) / ((1 - a) * (1 - 6 * u) + 2 * u)

        J0_exp = ((1 - u)^2) * (a + (1 - 2 * a) * (1 - u)^2) / (2 * N * G)
        J1_exp = ((1 - u)^2) * (1 - a) / (2 * N * G)

        J0_exp2 = 1 / (8 * N * u * c + 1)
        J1_exp2 = (1 - a) * (1 - 2 * u) / (8 * N * u * (1 - a + u) + (1 - a) * (1 - 6 * u) + 2 * u)

        lambda_1 = a*(1-1/(4*N)) + sqrt((a^2)*(1-1/(4*N))^2 + (1-2*a)*(1-1/(2*N)))
        lambda_2 = a*(1-1/(4*N)) - sqrt((a^2)*(1-1/(4*N))^2 + (1-2*a)*(1-1/(2*N)))
        
        vector_11 = 1-a
        vector_12 = -a*(1-1/(2*N)) + lambda_1
        
        vector_21 = 1-a
        vector_22 = -a*(1-1/(2*N)) + lambda_2
        
        eigen_vector_1 = [vector_11, vector_12]
        eigen_vector_2 = [vector_21, vector_22]
        
        M = [eigen_vector_1']

        exp_values_g[m] = (J0_exp, J1_exp)
        speed_g[m] = lambda_1 * ((1- u)^2)

        eq = zeros(Float64, 1)
        slope = zeros(Float64, 1)
        combinations = [(0.0, 0.0)]

        for (idx, (J0_val, J1_val)) in enumerate(combinations)
            J1_vec = fill(NaN, nGen)
            J0_vec = fill(NaN, nGen)

            J1_vec[1] = J1_val
            J0_vec[1] = J0_val

            for t in 2:nGen # iterate recursions, stop when reaching equilibrium value 
                J0_vec[t] = ((1-u)^2) * ((((1-m)^2) + m^2) * (((1/(2*N)) + ((1-1/(2*N)) * J0_vec[t-1]))) + (1 - ((1-m)^2) - m^2) * J1_vec[t-1])
                J1_vec[t] = ((1-u)^2) * (2*m*(1-m) * ((1/(2*N)) + (1-1/(2*N)) * J0_vec[t-1]) + (1 - 2*m*(1-m)) * J1_vec[t-1])

                if J0_vec[1] < J0_exp && J1_vec[1] < J1_exp
                    if J0_vec[t] >= J0_exp - 0.00001 && J1_vec[t] >= J1_exp - 0.00001 && eq[idx] == 0
                        eq[idx] = t
                    end
                elseif J0_vec[1] < J0_exp && J1_vec[1] > J1_exp
                    if J0_vec[t] >= J0_exp - 0.00001 && J1_vec[t] <= J1_exp + 0.00001 && eq[idx] == 0
                        eq[idx] = t
                    end
                elseif J0_vec[1] > J0_exp && J1_vec[1] < J1_exp
                    if J0_vec[t] <= J0_exp + 0.00001 && J1_vec[t] >= J1_exp - 0.00001 && eq[idx] == 0
                        eq[idx] = t
                    end
                elseif J0_vec[1] > J0_exp && J1_vec[1] > J1_exp
                    if J0_vec[t] <= J0_exp + 0.00001 && J1_vec[t] <= J1_exp + 0.00001 && eq[idx] == 0
                        eq[idx] = t
                    end
                end
                if eq[idx] == t
                    break
                end
            end

            EQ = Int(eq[idx]) # Time near-equilibrium was reached
            slope[idx] = abs(mean(diff(J1_vec[EQ-1:EQ]) ./ diff(J0_vec[EQ-1:EQ])))
            println(EQ)
            println("slope", slope)
            results_theoretical[(m)] = (
                J0_vec = J0_vec,
                J1_vec = J1_vec,
                EQ = EQ,
                slope = slope[idx])
        end

        # Calculate the mean slope for the current `m`
        slopeG = mean(skipmissing(slope))
        int = J1_exp - (slopeG * J0_exp) # Find the intercept from numerical data
    end
end


# CULTURAL EXPECTATION
results_culture_theoretical = Dict()
exp_values_c = Dict() 
speed_c = Dict()
begin
    for m in m_values
        a = ((1 - m)^2) + m^2
        G = 1 - (a * (1 - nu)^2) * (2 - 1/(N)) - (1 - 2 * a) * ((1 - nu)^4) * (1 - 1/(N))

        J0_exp_c = ((1 - nu)^2) * (a + (1 - 2 * a) * (1 - nu)^2) / (N * G)
        J1_exp_c = ((1 - nu)^2) * (1 - a) / (N * G)

        lambda_1_c = ((2 * N * a) - a + (4 * N^2 * a^2 - 8 * N^2 * a + 4 * N^2 - 4 * N * a^2 + 8 * N * a - 4 * N + a^2)^(1/2)) / (2*N)
        lambda_2_c = - (a - (2 * N * a) + (4 * N^2 * a^2 - 8 * N^2 * a + 4 * N^2 - 4 * N * a^2 + 8 * N * a - 4 * N + a^2)^(1/2)) / (2*N)

        vector_11 = 1-a
        vector_12 = -a*(1-1/(N))+lambda_1_c
        
        vector_21 = 1-a
        vector_22 = -a*(1-1/(N))+lambda_2_c

        eigen_vector_1_c = [vector_11, vector_12]
        eigen_vector_2_c = [vector_21, vector_22]
        
        M_c = [eigen_vector_1_c']

        exp_values_c[m] = (J0_exp_c, J1_exp_c)
        speed_c[m] = lambda_1_c * ((1-nu)^2)

        eq_c = zeros(Float64, 1)
        slope_c = zeros(Float64, 1)
        combinations = [(0.0, 0.0)]

        for (idx, (J0_val, J1_val)) in enumerate(combinations)
            J1_vec = fill(NaN, nGen)
            J0_vec = fill(NaN, nGen)
            J1_vec[1] = J1_val
            J0_vec[1] = J0_val

            for t in 2:nGen
                J0_vec[t] = ((1-nu)^2) * ((((1-m)^2) + m^2) * ((1/N) + ((1-(1/N)) * J0_vec[t-1])) + (1 - ((1-m)^2) - m^2) * J1_vec[t-1])
                J1_vec[t] = ((1-nu)^2) * (2*m*(1-m) * ((1/N) + (1-(1/N)) * J0_vec[t-1]) + (1 - 2*m*(1-m)) * J1_vec[t-1])

                if J0_vec[1] < J0_exp_c && J1_vec[1] < J1_exp_c
                    if J0_vec[t] >= J0_exp_c - 0.00001 && J1_vec[t] >= J1_exp_c - 0.00001 && eq_c[idx] == 0
                        eq_c[idx] = t
                    end
                elseif J0_vec[1] < J0_exp_c && J1_vec[1] > J1_exp_c
                    if J0_vec[t] >= J0_exp_c - 0.00001 && J1_vec[t] <= J1_exp_c + 0.00001 && eq_c[idx] == 0
                        eq_c[idx] = t
                    end
                elseif J0_vec[1] > J0_exp_c && J1_vec[1] < J1_exp_c
                    if J0_vec[t] <= J0_exp_c + 0.00001 && J1_vec[t] >= J1_exp_c - 0.00001 && eq_c[idx] == 0
                        eq_c[idx] = t
                    end
                elseif J0_vec[1] > J0_exp_c && J1_vec[1] > J1_exp_c
                    if J0_vec[t] <= J0_exp_c + 0.00001 && J1_vec[t] <= J1_exp_c + 0.00001 && eq_c[idx] == 0
                        eq_c[idx] = t
                    end
                end

                if eq_c[idx] == t
                    break
                end
            end

            EQ = Int(eq_c[idx]) 
            println(EQ)
            slope_c[idx] = abs(mean(diff(J1_vec[EQ-10:EQ]) ./ diff(J0_vec[EQ-10:EQ])))

            results_culture_theoretical[(m)] = (
                J0_vec = J0_vec,
                J1_vec = J1_vec,
                EQ = EQ,
                slope_c = slope_c[idx])
        end
        slopeC = mean(skipmissing(slope_c))
        intC = J1_exp_c - (slopeC * J0_exp_c) 
    end
end

# color palette, modified original :Spectral palette to avoid bright yellow 
my_colors = [RGBA(0.62, 0.004, 0.259, 1.0),
            RGBA(0.7275, 0.1235, 0.2845, 1.0),
            RGBA(0.835, 0.243, 0.31, 1.0),
            RGBA(0.8959999999999999, 0.33499999999999996, 0.2865, 1.0),
            RGBA(0.957, 0.427, 0.263, 1.0),
            RGBA(0.9744999999999999, 0.5545, 0.3215, 1.0),
            RGBA(0.992, 0.682, 0.38, 1.0),
            RGBA(0.994, 0.78, 0.4625, 1.0),
            RGBA(0.996, 0.878, 0.545, 1.0),
            RGBA(0.998, 0.939, 0.647, 1.0),  
            RGBA(0.88, 0.93, 0.62, 1.0), 
            RGBA(0.85, 0.91, 0.58, 1.0),
            RGBA(0.902, 0.961, 0.596, 1.0),
            RGBA(0.7865, 0.9139999999999999, 0.6194999999999999, 1.0), 
            RGBA(0.671, 0.867, 0.643, 1.0),
            RGBA(0.5355000000000001, 0.8140000000000001, 0.645, 1.0),
            RGBA(0.4, 0.761, 0.647, 1.0),
            RGBA(0.29800000000000004, 0.647, 0.694, 1.0),
            RGBA(0.196, 0.533, 0.741, 1.0),
            RGBA(0.2825, 0.4215, 0.688, 1.0),
            RGBA(0.369, 0.31, 0.635, 1.0),
            RGBA(133/255,60/255,120/255,1.0)]




# equilibrium plot :) 
begin
    color_palette = cgrad(my_colors)
    line_styles = [:solid, :dash, :dot, :dashdot, :longdash]
    p = plot(dpi=1000)
    u_ = round(u, digits = 8)
    nu_ = round(nu, digits = 4)
    for i in 1:length(m_values) 
        m_ = round(m_values[i], digits = 4)
        results_ = results[i][end-49999:end, :] # take last 50,000 time steps 
        J0_values= [mean(results_.J0Con1P1), mean(results_.J0Con2P1), mean(results_.J0Con3P1), mean(results_.J0Con4P1), mean(results_.J0Con5P1),
                    mean(results_.J0Con6P1), mean(results_.J0Con7P1), mean(results_.J0Con8P1),                    
                    mean(results_.J0Unb1), 
                    mean(results_.J0Con9P1), mean(results_.J0Con10P1), 
                    mean(results_.J0Con11P1), mean(results_.J0Con12P1), mean(results_.J0Con13P1), mean(results_.J0Con14P1), mean(results_.J0Con15P1),
                    mean(results_.J0Con16P1), mean(results_.J0Con17P1), mean(results_.J0Con18P1), mean(results_.J0Con19P1), mean(results_.J0Con20P1)]

        J1_values= [mean(results_.J1Con1), mean(results_.J1Con2), mean(results_.J1Con3), mean(results_.J1Con4), mean(results_.J1Con5),
                    mean(results_.J1Con6), mean(results_.J1Con7), mean(results_.J1Con8), 
                    mean(results_.J1Unb), 
                    mean(results_.J1Con9), mean(results_.J1Con10),
                    mean(results_.J1Con11), mean(results_.J1Con12), mean(results_.J1Con13), mean(results_.J1Con14), mean(results_.J1Con15),
                    mean(results_.J1Con16), mean(results_.J1Con17), mean(results_.J1Con18), mean(results_.J1Con19), mean(results_.J1Con20)]
        
        base_line_palette = cgrad(my_colors, length(J0_values), categorical = true) 

        for j in 1:20 # plot lines
           segment_color = base_line_palette[j]
           p = plot!([J0_values[j], J0_values[j+1]], [J1_values[j], J1_values[j+1]], linestyle = line_styles[i], linewidth = 2, alpha = 1,color=segment_color, legend = false, background_color=:transparent)
        end

        # plot scatter dots
        p = scatter!([mean(results_.J0Con1P1)], [mean(results_.J1Con1)], label="", color=color_palette[1], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black,
        legend=:outerright, fg_legend=:transparent, xlims = (0,1.1), ylims = (-0.1,1.1), aspect_ratio=:equal)        
        p = scatter!([mean(results_.J0Con2P1)], [mean(results_.J1Con2)], label="", color=color_palette[2], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con3P1)], [mean(results_.J1Con3)], label="", color=color_palette[3], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con4P1)], [mean(results_.J1Con4)], label="", color=color_palette[4], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con5P1)], [mean(results_.J1Con5)], label="", color=color_palette[5], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con6P1)], [mean(results_.J1Con6)], label="", color=color_palette[6], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con7P1)], [mean(results_.J1Con7)], label="", color=color_palette[7], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con8P1)], [mean(results_.J1Con8)], label="", color=color_palette[8], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)

        p = scatter!([mean(results_.J0Con9P1)], [mean(results_.J1Con9)], label="", color=color_palette[9], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con10P1)], [mean(results_.J1Con10)], label="", color=color_palette[10], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con11P1)], [mean(results_.J1Con11)], label="", color=color_palette[11], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con12P1)], [mean(results_.J1Con12)], label="", color=color_palette[12], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con13P1)], [mean(results_.J1Con13)], label="", color=color_palette[13], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con14P1)], [mean(results_.J1Con14)], label="", color=color_palette[14], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con15P1)], [mean(results_.J1Con15)], label="", color=color_palette[15], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con16P1)], [mean(results_.J1Con16)], label="", color=color_palette[16], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con17P1)], [mean(results_.J1Con17)], label="", color=color_palette[17], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con18P1)], [mean(results_.J1Con18)], label="", color=color_palette[18], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con19P1)], [mean(results_.J1Con19)], label="", color=color_palette[19], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black)
        p = scatter!([mean(results_.J0Con20P1)], [mean(results_.J1Con20)], label="", color=color_palette[20], markerstrokewidth = 1, markersize = 2.5, markerstrokecolor = :black, legend = false)

        p = xlabel!("\$Similarity \\ within  \\ a \\ population \\ J_0\$", xlabelfontsize = 10 )
        p = ylabel!("\$Similarity \\ between \\ populations \\ J_1\$", ylabelfontsize = 10)
    end  
    
    for i in 1:length(m_values) # plot equilibrium values of simulations and expectations
        m_ = round(m_values[i], digits = 4)
        results_ = results[i][end-49999:end, :] # take last 50,000 time steps 
        p = scatter!( [exp_values_c[m_values[i]][1]], [exp_values_c[m_values[i]][2]], color =:thistle, markerstrokewidth = 0, markersize = 8, marker =:star, legend = false)
        p = scatter!([mean(results_.J0Unb1)],[mean(results_.J1Unb)], color =:plum4, markerstrokewidth = 0, markersize = 4, marker =:diamond, legend = false)

        p = scatter!( [exp_values_g[m_values[i]][1]], [exp_values_g[m_values[i]][2]], color=:pink, markerstrokewidth = 0, markersize = 8, marker =:star, legend = false)
        p = scatter!([mean(results_.J0Gen1)], [mean(results_.J1Gen)], color=:palevioletred3, markerstrokewidth = 0, markersize = 4, marker =:diamond, legend = false)
    end

    results_ = results[5][end-49999:end, :] # take last 50,000 time steps 
    J0_values = [mean(results_.J0Con1P1), mean(results_.J0Con2P1), mean(results_.J0Con3P1),
                mean(results_.J0Con4P1), mean(results_.J0Con5P1), mean(results_.J0Con6P1),
                mean(results_.J0Con7P1), mean(results_.J0Con8P1), mean(results_.J0Unb1),
                mean(results_.J0Con9P1), mean(results_.J0Con10P1), mean(results_.J0Con11P1),
                mean(results_.J0Con12P1), mean(results_.J0Con13P1), mean(results_.J0Con14P1),
                mean(results_.J0Con15P1), mean(results_.J0Con16P1), mean(results_.J0Con17P1),
                mean(results_.J0Con18P1), mean(results_.J0Con19P1), mean(results_.J0Con20P1)]
    for (i, x) in enumerate(J0_values) # plot rug 
        plot!(p, [x, x], [-0.05, -0.06], linewidth = 2, color = color_palette[i])
    end
    x_star = mean(results_.J0Unb1)
    scatter!([x_star], [-0.05], marker = (:diamond, 3.8), color =:plum4, markerstrokewidth = 0) 
    x_star = mean(results_.J0Gen1)
    scatter!([x_star], [-0.05], marker = (:diamond, 3.8), color =:palevioletred3, markerstrokewidth = 0) 

    # save figures :) 
    savefig(p, "$dir/phaseplot_N_$(N)_m_$(m_values)_u_$(u_)_nu_($nu)_with_rug.png") 
    savefig(p, "$dir/phaseplot_N_$(N)_m_$(m_values)_u_$(u_)_nu_($nu)_with_rug.svg")     
end

begin # plot legend 
    line_styles = [:solid, :dash, :dot, :dashdot, :longdash]
    p_legend = plot(xlim=(0,0.1), ylim=(0,0.1), framestyle=:none, dpi = 200, grid=false,background_color=:transparent, fg_legend=:transparent)
    dummy_x = [0.1, 0.1]
    dummy_y = [0.1, 0.1]
    for i in 1:length(m_values)
        plot!(p_legend, dummy_x, dummy_y, linestyle=line_styles[i], color=:black, label="m=$(round(m_values[i], digits=4))", legendfontsize=14)
    end
    savefig(p_legend, "$dir/legend_migration_linestyles.png")
    savefig(p_legend, "$dir/legend_migration_linestyles.svg")
end