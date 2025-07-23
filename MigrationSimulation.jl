using Distributed
addprocs(200) # add worker processes for parallel computing
@everywhere using Random
@everywhere using Dates
@everywhere using StatsBase
@everywhere using Distributions
@everywhere using DelimitedFiles
@everywhere using DataFrames
@everywhere using CSV
@everywhere using SharedArrays

# FUNCTIONS # 
begin 
    # create two populations ; change for other initial J0 and J1 
    @everywhere function initialisePopulations!(data1, data2, N, gen1, gen2, unb, conf1, conf2, conf3, conf4, conf5, conf6, conf7, conf8, conf9, conf10, conf11, conf12, conf13, conf14, conf15, conf16, conf17, conf18, conf19, conf20)  
        data1[gen1,:] .=  collect(0.0001:0.0001:0.0001*N)
        data1[gen2,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data1[unb,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf1,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf2,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf3,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf4,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf5,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf6,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf7,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf8,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf9,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf10,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf11,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf12,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf13,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf14,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf15,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf16,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf17,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf18,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf19,:] .= collect(0.0001:0.0001:0.0001*N)
        data1[conf20,:] .= collect(0.0001:0.0001:0.0001*N)

        data2[gen1,:] .= collect(0.2001:0.0001:0.0001*(3*N))
        data2[gen2,:] .= collect(0.3001:0.0001:0.0001*(4*N))
        data2[unb,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf1,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf2,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf3,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf4,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf5,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf6,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf7,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf8,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf9,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf10,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf11,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf12,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf13,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf14,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf15,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf16,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf17,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf18,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf19,:] .= collect(0.1001:0.0001:0.0001*(2N))
        data2[conf20,:] .= collect(0.1001:0.0001:0.0001*(2N))
    end
    # as the name suggests, it initializes output containers
    @everywhere function initialize_output_containers(tsamps, calculateIBD, calculateVars, calculateFST) 
        J0gen1 = J0gen2 = J1gen = J0unb1 = J0unb2 = J1unb = nothing
        J0conf1p1 = J0conf1p2 = J1conf1 = nothing
        J0conf2p1 = J0conf2p2 = J1conf2 = nothing
        J0conf3p1 = J0conf3p2 = J1conf3 = nothing
        J0conf4p1 = J0conf4p2 = J1conf4 = nothing
        J0conf5p1 = J0conf5p2 = J1conf5 = nothing
        J0conf6p1 = J0conf6p2 = J1conf6 = nothing
        J0conf7p1 = J0conf7p2 = J1conf7 = nothing
        J0conf8p1 = J0conf8p2 = J1conf8 = nothing
        J0conf9p1 = J0conf9p2 = J1conf9 = nothing
        J0conf10p1 = J0conf10p2 = J1conf10 = nothing
        J0conf11p1 = J0conf11p2 = J1conf11 = nothing
        J0conf12p1 = J0conf12p2 = J1conf12 = nothing
        J0conf13p1 = J0conf13p2 = J1conf13 = nothing
        J0conf14p1 = J0conf14p2 = J1conf14 = nothing
        J0conf15p1 = J0conf15p2 = J1conf15 = nothing
        J0conf16p1 = J0conf16p2 = J1conf16 = nothing
        J0conf17p1 = J0conf17p2 = J1conf17 = nothing
        J0conf18p1 = J0conf18p2 = J1conf18 = nothing
        J0conf19p1 = J0conf19p2 = J1conf19 = nothing
        J0conf20p1 = J0conf20p2 = J1conf20 = nothing
        if calculateIBD
            J0gen1, J0gen2, J1gen = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0unb1, J0unb2, J1unb = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            
            J0conf1p1, J0conf1p2, J1conf1 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf2p1, J0conf2p2, J1conf2 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf3p1, J0conf3p2, J1conf3 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf4p1, J0conf4p2, J1conf4 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf5p1, J0conf5p2, J1conf5 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf6p1, J0conf6p2, J1conf6 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf7p1, J0conf7p2, J1conf7 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf8p1, J0conf8p2, J1conf8 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf9p1, J0conf9p2, J1conf9 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf10p1, J0conf10p2, J1conf10 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf11p1, J0conf11p2, J1conf11 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf12p1, J0conf12p2, J1conf12 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf13p1, J0conf13p2, J1conf13 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf14p1, J0conf14p2, J1conf14 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf15p1, J0conf15p2, J1conf15 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf16p1, J0conf16p2, J1conf16 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf17p1, J0conf17p2, J1conf17 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf18p1, J0conf18p2, J1conf18 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf19p1, J0conf19p2, J1conf19 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
            J0conf20p1, J0conf20p2, J1conf20 = zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps)), zeros(Float64, 1, length(tsamps))
        end

        VarSumTraits = nothing
        VarSumGen = nothing
        if calculateVars
            VarSumTraits = [Vector{Tuple{Int, Int, Int, Float64, Float64, Float64, Float64, Int}}(undef, length(tsamps)) for _ in 1:21]
            VarSumGen = Vector{Tuple{Int, Int, Int, Float64, Float64, Float64, Float64, Int}}(undef, length(tsamps))
        end

        ValsFST = nothing
        if calculateFST
            ValsFST = [Vector{Float64}(undef, length(tsamps)) for _ in 1:11] # modify this thing for fst 
        end
    
        return J0gen1, J0gen2, J1gen, J0unb1, J0unb2, J1unb,
               J0conf1p1, J0conf1p2, J1conf1, J0conf2p1, J0conf2p2, J1conf2,
               J0conf3p1, J0conf3p2, J1conf3, J0conf4p1, J0conf4p2, J1conf4,
               J0conf5p1, J0conf5p2, J1conf5, J0conf6p1, J0conf6p2, J1conf6,
               J0conf7p1, J0conf7p2, J1conf7, J0conf8p1, J0conf8p2, J1conf8,
               J0conf9p1, J0conf9p2, J1conf9, J0conf10p1, J0conf10p2, J1conf10, 
               J0conf11p1, J0conf11p2, J1conf11, J0conf12p1, J0conf12p2, J1conf12,
               J0conf13p1, J0conf13p2, J1conf13, J0conf14p1, J0conf14p2, J1conf14,
               J0conf15p1, J0conf15p2, J1conf15, J0conf16p1, J0conf16p2, J1conf16,
               J0conf17p1, J0conf17p2, J1conf17, J0conf18p1, J0conf18p2, J1conf18,
               J0conf19p1, J0conf19p2, J1conf19, J0conf20p1, J0conf20p2, J1conf20,
               VarSumTraits, VarSumGen, ValsFST    
    end
    # reproduction: first makes a genetic pool, them samples from that pool     
    @everywhere function reproduction!(data, N)
        genome = vcat(data[1,:], data[2,:])
        parents = rand(1:2 * N, 2 * N)
        sampled_alleles = genome[parents]
        data[1, :] = sampled_alleles[1:N]
        data[2, :] = sampled_alleles[N + 1:2 *N]
    end
    # unbiased transmission
    @everywhere function unbiased(culpool, N)
            indices = rand(1:N, N)
            return culpool[indices]
    end 
    # frequency dependent transmission, although I named the function conformity it in fact works for both anticonformity or conformity :) 
    @everywhere function conformity(culpool, θ, N)
            uniqueTraits = sort(unique(culpool))
            culFreq = [count(x -> x == trait, culpool) / length(culpool) for trait in uniqueTraits] 
            culFreq = (culFreq .^ (1 + θ)) / sum(culFreq .^ (1 + θ))
            culconf = sample(uniqueTraits, Weights(culFreq), N)
            return culconf
    end 
    # same as conformity but iterates over many traits to not have to apply the function to each trait
    @everywhere function conformityManyTraits!(data, freqdep_indices, θs, N)
            for i in 1:length(freqdep_indices)
                culpool = data[freqdep_indices[i], :]
                uniqueTraits = sort(unique(culpool))
                culFreq = [count(x -> x == trait, culpool) / length(culpool) for trait in uniqueTraits]
                culFreq = (culFreq .^ (1 + θs[i])) / sum(culFreq .^ (1 + θs[i]))
                culconf = sample(uniqueTraits, Weights(culFreq), N)
                data[freqdep_indices[i], :] .= culconf
            end
    end 
    # implements mutations 
    @everywhere function mutate!(data, mutation_rate, population_size) 
                mutnum = rand(Binomial(2 * population_size, mutation_rate)) # number
                if mutnum > 0
                    mutants = rand(1:population_size, mutnum) # indices
                    for i in 1:mutnum
                       selected_gen = rand(Bool) ? :1 : :2
                       data[selected_gen,mutants[i]] = rand(mutnum)[i]
                    end
                end
    end
    # implements innovations 
    @everywhere function innovate!(data, innovation_rate, population_size,row) 
                innum = rand(Binomial(population_size, innovation_rate)) # number
                if innum > 0
                    innovators = rand(1:population_size, innum) # indices
                    for i in 1:innum
                        data[row, innovators[i]] = rand(innum)[i]
                    end
                end
    end    
    # implements migration 
    @everywhere function migrate!(data1, data2, migration_rate, population_size)
                mignum = rand(Binomial(population_size, migration_rate))
                if mignum > 0
                    migrants1 = rand(1:population_size, mignum)
                    migrants2 = rand(1:population_size, mignum)
                    moving1 = copy(data1[:,migrants1])
                    moving2 = copy(data2[:,migrants2])
                    data1[:,migrants1] .= moving2
                    data2[:,migrants2] .= moving1
                end
    end
    # calculates J0 and J1 for genetic trait 
    @everywhere function calculateGenetics!(data1, data2, ts, J0gen1, J0gen2, J1gen, N)
        testa = rand(1: 2 * N, 2 * N)
        testb = rand(1: 2 * N, 2 * N)
        testc = rand(1: 2 * N, 2 * N)
        testd = rand(1: 2 * N, 2 * N)
        genome1 = vcat(data1[1,:], data1[2,:])
        genome2 = vcat(data2[1,:], data2[2,:])
        J0gen1[1, ts] = sum(genome1[testa] .== genome1[testb]) / (2 * N)
        J0gen2[1, ts] = sum(genome2[testa] .== genome2[testb]) / (2 * N)
        J1gen[1, ts] = sum(genome1[testc] .== genome2[testd]) / (2 * N)
    end
    # calculates J0 and J1 for cultural traits 
    @everywhere function calculateCulture!(data1, data2, ts, J0_1, J0_2, J1, trait, N)
            testa = rand(1: N, N) 
            testb = rand(1: N, N) 
            testc = rand(1: N, N) 
            testd = rand(1: N, N) 
            J0_1[1, ts] = sum(data1[trait,:][testa] .== data1[trait,:][testb]) / N 
            J0_2[1, ts] = sum(data2[trait,:][testa] .== data2[trait,:][testb]) / N
            J1[1, ts] = sum(data1[trait,:][testc] .== data2[trait,:][testd]) / N
    end
   # calculate for cultural traits number of variants, shared variants, most abundant trait and total number of traits in the meta population  
    @everywhere function calculateCombinedSummary(data1, data2, trait)
            trait = trait + 2 # move it two indices (skip genetic trait, two because it is diploid)
            combined_repertoire = vcat(data1[trait, :], data2[trait, :])
            uniqueV = unique(combined_repertoire)
            
            variantC_P1 = [count(data1[trait, :] .== v) for v in uniqueV]
            variantC_P2 = [count(data2[trait, :] .== v) for v in uniqueV]
            
            num_variants_data1 = length(unique(data1[trait, :]))
            num_variants_data2 = length(unique(data2[trait, :]))
            
            shared_variants = length(intersect(unique(data1[trait, :]), unique(data2[trait, :])))
            
            most_abundant_variant_data1 = uniqueV[argmax(variantC_P1)]
            most_abundant_count_data1 = maximum(variantC_P1)/N
            
            most_abundant_variant_data2 = uniqueV[argmax(variantC_P2)]
            most_abundant_count_data2 = maximum(variantC_P2)/N
            
            #same_most_abundant_variant = Int(most_abundant_variant_data1 == most_abundant_variant_data2)

            return  num_variants_data1, num_variants_data2, shared_variants, 
                    most_abundant_variant_data1, most_abundant_count_data1, 
                    most_abundant_variant_data2, most_abundant_count_data2, length(uniqueV)
    end

    # calculate for genetic trait number of variants, shared variants, most abundant trait and total number of traits in the meta population  
    @everywhere function calculateCombinedSummaryGens(data1, data2)
            genome1 = vcat(data1[1,:], data1[2,:])
            genome2 = vcat(data2[1,:], data2[2,:])
            combined_repertoire = vcat(genome1, genome2)
            uniqueV = unique(combined_repertoire)
            
            variantC_P1 = [count(genome1 .== v) for v in uniqueV]
            variantC_P2 = [count(genome2 .== v) for v in uniqueV]
            
            num_variants_data1 = length(unique(genome1))
            num_variants_data2 = length(unique(genome2))
            
            shared_variants = length(intersect(unique(genome1), unique(genome2)))
            
            most_abundant_variant_data1 = uniqueV[argmax(variantC_P1)]
            most_abundant_count_data1 = maximum(variantC_P1)/(2*N)
            
            most_abundant_variant_data2 = uniqueV[argmax(variantC_P2)]
            most_abundant_count_data2 = maximum(variantC_P2)/(2*N)

            #same_most_abundant_variant = Int(most_abundant_variant_data1 == most_abundant_variant_data2)

            return  num_variants_data1, num_variants_data2, shared_variants, 
                    most_abundant_variant_data1, most_abundant_count_data1, 
                    most_abundant_variant_data2, most_abundant_count_data2, length(uniqueV) 
    end

    # Fst calculation implemented as in Mesoudi 2018 Supplementary S1 Methods 
    @everywhere function calculateFSTstat(data1, data2, trait) 
        index = trait + 2 
        # var total
        combined_repertoire = vcat(data1[index, :], data2[index, :])
        uniqueV = unique(combined_repertoire)
        total_individuals  = size(data1, 2) + size(data2, 2)
        mean_freqs = [(count(==(v), data1[index, :]) + count(==(v), data2[index, :])) / total_individuals for v in uniqueV ]
        var_total = 1 - sum(mean_freqs .^ 2)

        # var within 
        unique1 = unique(data1[index, :])
        frequencies1 = [count(==(v), data1[index, :]) / size(data1, 2) for v in unique(data1[index, :])]
        var_within1 = 1 - sum(frequencies1 .^2)


        unique2 = unique(data2[index, :])
        frequencies2= [count(==(v), data2[index, :]) / size(data2, 2) for v in unique(data2[index, :])]
        var_within2 = 1 - sum(frequencies2 .^2)        

        var_within = (var_within1 + var_within2) / 2

        fst = (var_total - var_within) / var_total
        return fst
    end

    # Replicate run 
    @everywhere function run_replicate!(r, N, u, nu, θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20, m, nGen, tsamps, numTraits)
        println("replicate $r m $m nu $nu")
        gen1 = 1    # allele 1
        gen2 = 2    # allele 2
        unb = 3     # cultural trait unbiased
        conf1 = 4    # cultural trait freq-dep
        conf2 = 5    # cultural trait freq-dep
        conf3 = 6    # cultural trait freq-dep
        conf4 = 7    # cultural trait freq-dep
        conf5 = 8    # cultural trait freq-dep
        conf6 = 9    # cultural trait freq-dep
        conf7 = 10    # cultural trait freq-dep
        conf8 = 11    # cultural trait freq-dep
        conf9 = 12    # cultural trait freq-dep
        conf10 = 13    # cultural trait freq-dep
        conf11 = 14    # cultural trait freq-dep
        conf12 = 15    # cultural trait freq-dep
        conf13 = 16    # cultural trait freq-dep
        conf14 = 17    # cultural trait freq-dep
        conf15 = 18    # cultural trait freq-dep
        conf16 = 19    # cultural trait freq-dep
        conf17 = 20    # cultural trait freq-dep
        conf18 = 21    # cultural trait freq-dep
        conf19 = 22    # cultural trait freq-dep
        conf20 = 23    # cultural trait freq-dep

        
        freqdep_indices = [conf1, conf2, conf3, conf4, conf5, conf6, conf7, conf8, conf9, conf10, conf11, conf12, conf13, conf14, conf15, conf16, conf17, conf18, conf19, conf20]
        conf_indices = [conf1, conf2, conf3, conf4, conf5, conf6, conf7, conf8, conf9, conf10, conf11, conf12, conf13, conf14, conf15, conf16, conf17, conf18, conf19, conf20]
        θs = [θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20]
        data1 = Array{Float64}(undef, numTraits + 1, N)
        data2 = Array{Float64}(undef, numTraits +1, N)

        begin # output containers 
            J0gen1, J0gen2, J1gen, J0unb1, J0unb2, J1unb,
            J0conf1p1, J0conf1p2, J1conf1, J0conf2p1, J0conf2p2, J1conf2,
            J0conf3p1, J0conf3p2, J1conf3, J0conf4p1, J0conf4p2, J1conf4,
            J0conf5p1, J0conf5p2, J1conf5, J0conf6p1, J0conf6p2, J1conf6,
            J0conf7p1, J0conf7p2, J1conf7, J0conf8p1, J0conf8p2, J1conf8, 
            J0conf9p1, J0conf9p2, J1conf9, J0conf10p1, J0conf10p2, J1conf10, 
            J0conf11p1, J0conf11p2, J1conf11, J0conf12p1, J0conf12p2, J1conf12,
            J0conf13p1, J0conf13p2, J1conf13, J0conf14p1, J0conf14p2, J1conf14,
            J0conf15p1, J0conf15p2, J1conf15, J0conf16p1, J0conf16p2, J1conf16,
            J0conf17p1, J0conf17p2, J1conf17, J0conf18p1, J0conf18p2, J1conf18,
            J0conf19p1, J0conf19p2, J1conf19, J0conf20p1, J0conf20p2, J1conf20,
            VarSumTraits, VarSumGen, ValsFST = initialize_output_containers(tsamps, calculateIBD, calculateVars, calculateFST)

            conf_params = [(J0conf1p1, J0conf1p2, J1conf1), (J0conf2p1, J0conf2p2, J1conf2), 
                            (J0conf3p1, J0conf3p2, J1conf3), (J0conf4p1, J0conf4p2, J1conf4), 
                            (J0conf5p1, J0conf5p2, J1conf5), (J0conf6p1, J0conf6p2, J1conf6),
                            (J0conf7p1, J0conf7p2, J1conf7), (J0conf8p1, J0conf8p2, J1conf8),
                            (J0conf9p1, J0conf9p2, J1conf9), (J0conf10p1, J0conf10p2, J1conf10),
                            (J0conf11p1, J0conf11p2, J1conf11), (J0conf12p1, J0conf12p2, J1conf12),
                            (J0conf13p1, J0conf13p2, J1conf13), (J0conf14p1, J0conf14p2, J1conf14),
                            (J0conf15p1, J0conf15p2, J1conf15), (J0conf16p1, J0conf16p2, J1conf16),
                            (J0conf17p1, J0conf17p2, J1conf17), (J0conf18p1, J0conf18p2, J1conf18),
                            (J0conf19p1, J0conf19p2, J1conf19), (J0conf20p1, J0conf20p2, J1conf20)]
        end
        # initialise populations default J0 = 0 ; J1 = 0 max heterogeneity
        initialisePopulations!(data1, data2, N, gen1, gen2, unb, conf1, conf2, conf3, conf4, conf5, conf6, conf7, conf8, conf9, conf10, conf11, conf12, conf13, conf14, conf15, conf16, conf17, conf18, conf19, conf20) #ant1, ant2, ant3, ant4, ant5)
        ts = 0 # tsampling 
        for t in 1:nGen
            # samples for output 
            if t in tsamps
                ts += 1
                # sample IBD
                if calculateIBD == true 
                    calculateGenetics!(data1, data2, ts, J0gen1, J0gen2, J1gen, N)
                    calculateCulture!(data1, data2, ts, J0unb1, J0unb2, J1unb, unb, N)
                    for (i, (param1, param2, J1param)) in enumerate(conf_params)
                        calculateCulture!(data1, data2, ts, param1, param2, J1param, conf_indices[i], N)
                    end
                end 
                # sample VAR SUMS 
                if calculateVars == true #&& r in 1:1000
                    VarSumGen[ts] =  calculateCombinedSummaryGens(data1, data2)
                    for i in 1:21
                        VarSumTraits[i][ts] = calculateCombinedSummary(data1, data2, i)
                    end
                end 

                 # sample FST
                if calculateFST == true #&& r in 1:1000
                    for i in 1:11 # modify this depending on the fst 
                        ValsFST[i][ts] = calculateFSTstat(data1, data2, i)
                    end
                end 
            end   
    
            # reproduction
            reproduction!(data1, N)
            reproduction!(data2, N)

            # unbiased transmission 
            data1[unb,:] = unbiased(data1[unb, :], N)
            data2[unb,:] = unbiased(data2[unb, :], N)

            # conformity
            conformityManyTraits!(data1, freqdep_indices, θs, N)
            conformityManyTraits!(data2, freqdep_indices, θs, N)

            # mutate
            mutate!(data1, u, N)
            mutate!(data2, u, N)

            # innovate
            innovate!(data1, nu, N, 3)
            innovate!(data2, nu, N, 3)

            for i in freqdep_indices
                innovate!(data1, nu, N, i)
                innovate!(data2, nu, N, i)
            end

            # migrate
            migrate!(data1, data2, m, N)

        end
 

        # transform output containers to write individual files depending on which output flags where set as true 
        IBDmatrix = nothing
        if calculateIBD
            IBDmatrix = hcat(J0gen1', J0gen2', J1gen', J0unb1', J0unb2', J1unb',
            J0conf1p1', J0conf1p2', J1conf1', J0conf2p1', J0conf2p2', J1conf2',
            J0conf3p1', J0conf3p2', J1conf3', J0conf4p1', J0conf4p2', J1conf4',
            J0conf5p1', J0conf5p2', J1conf5', J0conf6p1', J0conf6p2', J1conf6', 
            J0conf7p1', J0conf7p2', J1conf7', J0conf8p1', J0conf8p2', J1conf8',
            J0conf9p1', J0conf9p2', J1conf9', J0conf10p1', J0conf10p2', J1conf10',
            J0conf11p1', J0conf11p2', J1conf11', J0conf12p1', J0conf12p2', J1conf12',
            J0conf13p1', J0conf13p2', J1conf13', J0conf14p1', J0conf14p2', J1conf14',
            J0conf15p1', J0conf15p2', J1conf15', J0conf16p1', J0conf16p2', J1conf16',
            J0conf17p1', J0conf17p2', J1conf17', J0conf18p1', J0conf18p2', J1conf18',
            J0conf19p1', J0conf19p2', J1conf19', J0conf20p1', J0conf20p2', J1conf20')
            if getIBDFiles 
                df = DataFrame(IBDmatrix, :auto)
                IBDFile = joinpath(dir, "J0_1_J1_0_replicate_$r.csv")
                CSV.write(IBDFile, df, header=false)
            end
        end

        VarsSummary = nothing
        if calculateVars 
            expColsCult = [reduce(hcat, [collect(v) for v in VarSumTrait]) for VarSumTrait in VarSumTraits]
            summaryVarsCult = transpose(vcat(expColsCult...))
            expColsGen = [reduce(hcat, [collect(v) for v in VarSumG]) for VarSumG in VarSumGen]
            summaryVarsGen = vcat(expColsGen...)
            VarsSummary = hcat(summaryVarsGen,summaryVarsCult)
            if getVarSumFiles 
                sum_VS = DataFrame(VarsSummary, :auto)
                VarSummaryFile = joinpath(dir, "VarSumAll_replicate_$r.csv")                
                CSV.write(VarSummaryFile, sum_VS, header=false)
            end
        end

        FSTValues = nothing
        if calculateFST
            FSTValues = hcat(ValsFST...)
            if getFSTFiles
                dfFST = DataFrame(ValsFST, :auto)
                FSTFile = joinpath(dir, "FST_replicate_$r.csv")    
                CSV.write(FSTFile, dfFST, header=false)
            end
        end

        return (IBDmatrix, VarsSummary, FSTValues)
    end


    @everywhere function run_simulation(reps, tsamps, N, u, nu, θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20, m, nGen, numTraits)
     
        IBDsum = getSumIBD ? SharedArray{Float64}(zeros((length(tsamps), 3*numTraits))) : nothing
        VarsSum = getSumVarSum ? SharedArray{Float64}(zeros((length(tsamps), 8*numTraits))) : nothing
        FSTSum = getSumFST ? SharedArray{Float64}(zeros((length(tsamps),11))) : nothing
        # run in parallel
        @sync begin
            @distributed for r in 1:reps
                replicate_result = run_replicate!(r, N, u, nu, θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20, m, nGen, tsamps, numTraits)

                if getSumIBD
                    IBDsum .+= replicate_result[1]
                end

                if getSumVarSum
                    VarsSum .+= replicate_result[2]
                end

                if getSumFST
                    FSTSum .+= replicate_result[3]
                end
            end
        end

        if getSumIBD 
            trait_names = [:J0Gen1, :J0Gen2, :J1Gen, :J0Unb1, :J0Unb2, :J1Unb, # modify here in case there are more or less traits evolving under conformity  
                           :J0Con1P1, :J0Con1P2, :J1Con1, :J0Con2P1, :J0Con2P2, 
                           :J1Con2, :J0Con3P1, :J0Con3P2, :J1Con3, :J0Con4P1, 
                           :J0Con4P2, :J1Con4, :J0Con5P1, :J0Con5P2, :J1Con5, 
                           :J0Con6P1, :J0Con6P2, :J1Con6, :J0Con7P1, :J0Con7P2, :J1Con7, 
                           :J0Con8P1, :J0Con8P2, :J1Con8, :J0Con9P1, :J0Con9P2, :J1Con9,
                           :J0Con10P1, :J0Con10P2, :J1Con10, :J0Con11P1, :J0Con11P2, :J1Con11,
                           :J0Con12P1, :J0Con12P2, :J1Con12, :J0Con13P1, :J0Con13P2, :J1Con13,
                           :J0Con14P1, :J0Con14P2, :J1Con14, :J0Con15P1, :J0Con15P2, :J1Con15,
                           :J0Con16P1, :J0Con16P2, :J1Con16, :J0Con17P1, :J0Con17P2, :J1Con17,
                           :J0Con18P1, :J0Con18P2, :J1Con18, :J0Con19P1, :J0Con19P2, :J1Con19,
                           :J0Con20P1, :J0Con20P2, :J1Con20]
            df = DataFrame(IBDsum, :auto)
            rename!(df, trait_names)
            CSV.write(joinpath(summary_dir, "sum_IBD_N_$(N)_m_$(m)_u_$(u)_nu_$(nu)_reps_$(reps)_thetaθ1$(θ1)toθ11$(θ20).csv"), df)
        end

        if getSumVarSum
            summary_names = [:nVarP1Gen, :nVarP2Gen, :shVarGen, :abVarP1Gen, :fabVarP1Gen, :abVarP2Gen, :fabVarP2Gen, :tVarGen, # likewise, modify here in case there are more or less traits evolving under conformity  
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
            df = DataFrame(VarsSum, :auto)
            rename!(df, summary_names)
            CSV.write(joinpath(summary_dir, "sum_VariantsSummary_N_$(N)_m_$(m)_u_$(u)_nu_$(nu)_reps_$(reps)_thetaθ1$(θ1)toθ11$(θ20).csv"),df) 
        end

        if getSumFST 
            df = DataFrame(FSTSum, :auto)
            CSV.write(joinpath(summary_dir, "sum_FST_N_$(N)_m_$(m)_u_$(u)_nu_$(nu)_reps_$(reps)_thetaθ1$(θ1)toθ11$(θ20).csv"), df, header=false)
        end


    end


    
end 

# Initial conditions 
@everywhere begin
        N = 10^3                                       # popualtion size
        u = 10^-4                                      # mutation rate
        nu = 0.0005                                    # innovation rate
        m = 0.002                                      # migration rate
        nGen = 500000                                  # number of generations
        reps = 1000                                    # number of replicates 
        tsamps = [nGen-49999:nGen ÷ nGen:nGen-1; nGen] # change here for when to sample 
        θ1 = 0.001                                     # θ values 
        θ2 = 0.002                                     # if you add or remove θ values, make sure you also modify the output containers for the correct size 
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
        # total number of traits (1 genetic (diploid), 1 unbiased, 20 frequency dependent)
        numTraits = 22
        # output flags 
        calculateIBD = false  # calculate similarity (J0 and J1)
        calculateVars = true  # calculate diveristy metrics (shared, total variants...)
        calculateFST = false  # calculate Fst
        getIBDFiles = false   # write ibd individual files 
        getSumIBD = false     # write the sum of the replicates to get mean trajectory 
        getVarSumFiles = true # write shared vars, most abundant vars ...
        getSumVarSum = true   # get sum of all replicate files of VarSum
        getFSTFiles = false   # get Fst for individual runs
        getSumFST = false     # get sum of Fst files
end 

# MAIN 
@everywhere timestamp = Dates.format(now(), "yy-mm-dd")
@everywhere dir = "./RESULTS/N_$(N)_m_$(m)_u_$(u)_nu_$(nu)_nGen_$(nGen)_reps_$(reps)_θsweep_$(timestamp)" # output directory for individual files 
mkdir(dir) # create output directory 
summary_dir = "./SUMMARY_FILES" # output directory for summary files

run_simulation(reps, tsamps, N, u, nu, θ1, θ2, θ3, θ4, θ5, θ6, θ7, θ8, θ9, θ10, θ11, θ12, θ13, θ14, θ15, θ16, θ17, θ18, θ19, θ20, m, nGen, numTraits) # run simulations


