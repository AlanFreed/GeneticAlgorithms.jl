module example1

# To use this code, you will need to download the genetic algorithm at:
# using Pkg
# Pkg.add(url = "https://github.com/AlanFreed/GeneticAlgorithms.jl")

using
    CairoMakie,       # Pixel based figure construction.
    GeneticAlgorithm,
    Statistics

export
    Ex1Species,
    runModel
#=
--------------------------------------------------------------------------------
=#

struct Ex1Species <: AbstractSpecies
    nExp::Int64                     number of experiments
    nCtl::Vector{Int64}             number of control  variables per experiment
    nRes::Vector{Int64}             number of response variables per experiment
    nPts::Vector{Int64}             number of datum points per experiment
    ctrl::Vector{Matrix{Float64}}   control  data for  each experiment
    resp::Vector{Matrix{Float64}}   response data from each experiment
    stdR::Vector{Vector{Float64}}   standard deviation in the responses per exp

    # constructor

    function Ex1Species()

        nExp = 1
        nCtl = Vector{Int64}(undef, nExp)
        nCtl[1] = 1     # i.e., strain
        nRes = Vector{Int64}(undef, nExp)
        nRes[1] = 1     # i.e., stress
        nPts = Vector{Int64}(undef, nExp)
        nPts[1] = 20
        ctrl = Vector{Matrix{Float64}}(undef, nExp)
        mtxC = Matrix{Float64}(undef, nCtl[1], nPts[1])
        mtxC[1,1]  = 0.1
        mtxC[1,2]  = 0.2
        mtxC[1,3]  = 0.3
        mtxC[1,4]  = 0.4
        mtxC[1,5]  = 0.5
        mtxC[1,6]  = 0.6
        mtxC[1,7]  = 0.7
        mtxC[1,8]  = 0.8
        mtxC[1,9]  = 0.9
        mtxC[1,10] = 1.0
        mtxC[1,11] = 1.1
        mtxC[1,12] = 1.2
        mtxC[1,13] = 1.3
        mtxC[1,14] = 1.4
        mtxC[1,15] = 1.5
        mtxC[1,16] = 1.6
        mtxC[1,17] = 1.7
        mtxC[1,18] = 1.8
        mtxC[1,19] = 1.9
        mtxC[1,20] = 2.0
        ctrl[1] = mtxC
        resp = Vector{Matrix{Float64}}(undef, nExp)
        mtxR = Matrix{Float64}(undef, nRes[1], nPts[1])
        mtxR[1,1]  = 0.050
        mtxR[1,2]  = 0.111
        mtxR[1,3]  = 0.193
        mtxR[1,4]  = 0.290
        mtxR[1,5]  = 0.349
        mtxR[1,6]  = 0.450
        mtxR[1,7]  = 0.559
        mtxR[1,8]  = 0.622
        mtxR[1,9]  = 0.744
        mtxR[1,10] = 0.835
        mtxR[1,11] = 1.032
        mtxR[1,12] = 1.144
        mtxR[1,13] = 1.266
        mtxR[1,14] = 1.396
        mtxR[1,15] = 1.409
        mtxR[1,16] = 1.494
        mtxR[1,17] = 1.625
        mtxR[1,18] = 1.675
        mtxR[1,19] = 1.700
        mtxR[1,20] = 1.710
        resp[1] = mtxR
        stdR = Vector{Vector{Float64}}(undef, nExp)
        stdV = Vector{Float64}(undef, nRes[1])
        stdV[1] = std(mtxR[1,:])
        stdR[1] = stdV
        new(nExp, nCtl, nRes, nPts, ctrl, resp, stdR)
    end
end # Ex1Species

function _ex1Model(ϵ::Float64, θ::Vector{Float64})::Float64
    # NOTE: the model's parameters are to index from 2,
    # as θ[1] is an internal parameter assigned by the optimizer.
    a = θ[2]
    b = θ[3]
    E = θ[4]
    σ = E * sinh(a*ϵ) / (a * cosh(a*ϵ) - b * sinh(a*\epsilon))
    return σ
end _ex1Model

function runModel(s::Ex1Species, θ::Vector{Float64})::Vector{Matrix{Float64}}
    # Retrieve the data held within MySpecies s.
    nExp = experiments(s)
    nRes = responsesPerExp(s)
    nPts = dataPointsPerExp(s)
    ctrl = controls(s)
    # This function returns the model's response as an array of dimension:
    #   [nExp] x [nRes x nPts]  indexed as  [i][j,k]
    #       nExp    number of experiments the model is to be fit against
    #       nRes    number of  responses  for each experiment  i
    #       nPts    number of data points for each experiment  i
    # Create the returned array that is to hold the model's response.
    modR = Vector{Matrix{Float64}}(undef, nExp)
    for exp in 1:nExp
        mtxR = Matrix{Float64}(undef, nRes[exp], nPts[exp])
        modR[exp] = mtxR
    end
    # Run your model subject to Ex1Species controls, viz., subject to ctrl.
    for i in 1:nExp
        for j in 1:nRes[i]
            for k in 1:nPts[i]
                ϵ = ctrl[i][j,k]
                σ = _ex1Model(ϵ, θ)
                modR[i][j,k] = σ
            end
        end
    end
    return modR
end
