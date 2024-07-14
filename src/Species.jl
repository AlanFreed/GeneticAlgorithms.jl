"""
This genetic algorithm internally calls the function:

function _evaluate(s::AbstractSpecies, c::Creature)
    ...
    θ = parameters(c)
    modelResponse = runModel(s, θ)
    ...
end

To use this genetic algorithm, a user must create their own data structure that
extends (or inherits) a base (or parent's) data structure described by:

@base struct AbstractSpecies
    nExp::Int                   # number of experiments
    nCtl::Vector{Int}           # number of control  variables per experiment
    nRes::Vector{Int}           # number of response variables per experiment
    nPts::Vector{Int}           # number of experimental data  per experiment
    ctrl::Vector{Matrix{Real}}  # control data:  [exp] x [nCtl[exp] x nPts[exp]]
    resp::Vector{Matrix{Real}}  # response data: [exp] x [nRes[exp] x nPts[exp]]
    stdR::Vector{Vector{Real}}  # standard deviations in the response variables
                                #    per experiment: [exp] x [nRes[exp]]
end

Arrays nCtl, nRes and nPts index as [i], where i ∈ [1, nExp].
Arrays ctrl and resp index as  [i][j,k], where i ∈ [1, nExp],
    j ∈ [1, nCtl[i] or nRes[i]] and k ∈ [1, nPts[i]].
Array stdR indexes as [i][j], where i ∈ [1, nExp] and j ∈ [1, nRes[i]].

Each extension (or child) type may introduce additional fields of data, if
required.  Typically, it won't, so a child with constructor might look like:

@extend struct MySpecies <: AbstractSpecies
    # internal constructor
    function MySpecies(nExp::Int, nCtl::Vector{Int}, nRes::Vector{Int}, nPts::Vector{Int}, ctrl::Vector{Matrix{Real}}, resp::Vector{Matrix{Real}})

        stdR = Vector{Vector{Real}}(undef, nExp)
        for exp = 1:nExp
            stdR[exp] = Vector{Real}(undef, nRes[exp])
            for res = 1:nRes[exp]
                stdR[exp][res] = std(resp[exp][res, :])
            end
        end

        new(nExp, nCtl, nRes, nPts, ctrl, resp, stdR)
    end
end

which 'inherits' fields: nExp, nCtl, nRes, nPts, ctrl, resp and stdR from its
parent type AbstractSpecies via the package ConcreteAbstractions.

All of these fields will be accessed by the genetic algorithm.

Along with the user's data structure, the user must supply their own function
(a template follows):

function runModel(s::MySpecies, θ::Vector{Real})::Vector{Matrix{Real}}
    # NOTE: the parameters (i.e., θ) for the model being fit start at index 2,
    # as θ[1] is an internal parameter assigned by the optimizer, viz.,
    # it is the 'p' for a p-norm of error.

    # This function returns the model's response as a vector of matrices:
    #   [exp] x [nRes[exp] x nPts[exp]]  indexed as  [i][j,k]  where
    #       exp         experiment number with values  exp ∈ [1, nExp]
    #       nRes[exp]   number of  responses  for experiment  exp
    #       nPts[exp]   number of data points for experiment  exp

    # Create the returned array that is to hold the model's response.
    modelResponse = Vector{Matrix{Real}}(undef, s.nExp)
    for exp in 1:s.nExp
        responseMatrix = Matrix{Real}(undef, s.nRes[exp], s.nPts[exp])
        modelResponse[exp] = responseMatrix
    end

    # Run your model subject to MySpecies s' controls, viz., subject to s.ctrl.
    # This model can be a function, a differential equation, an integral
    # equation, or whatever else the user may choose to implement here.  Be
    # aware that here is where the greatest expense of this algorithm resides.
    # ...
    for i in 1:s.nExp
        # ...
        for j in 1:s.nRes[i]
            # ...
            for k in 1:s.nPts[i]
                # ...
            end
            # ...
        end
        # ...
    end
    # ...
    return modelResponse
end
"""
@base struct AbstractSpecies
    nExp::Int                   # number of experiments
    nCtl::Vector{Int}           # number of control  variables per experiment
    nRes::Vector{Int}           # number of response variables per experiment
    nPts::Vector{Int}           # number of experimental data  per experiment
    ctrl::Vector{Matrix{Real}}  # control data:  [exp] x [nCtl[exp] x nPts[exp]]
    resp::Vector{Matrix{Real}}  # response data: [exp] x [nRes[exp] x nPts[exp]]
    stdR::Vector{Vector{Real}}  # standard deviations of the response variables
                                #    per experiment: [exp] x [nRes[exp]]
end

# Define a runModel function, using multiple dispatch to select the right one.

function runModel(s::AbstractSpecies, θ::Vector{Real})::Vector{Matrix{Real}}
    modelResponse = Vector{Matrix{Real}}(undef, 0)
    return modelResponse
end

#=
ConcreteAbstractions.md reads (edited for clarity):

Install with:

Pkg.add(url = "https://github.com/tbreloff/ConcreteAbstractions.jl")
using ConcreteAbstractions

The abstract definition:

@base struct AbstractFoo{T}
    a
    b::Int
    c::T
    d::Vector{T}
end

is really more like:

abstract struct AbstractFoo end
ConcreteAbstractions._base_types[:AbstractFoo] = ([:T], :(begin; a; b::Int; c::T; d::Vector{T}; end))

and the child definition:

@extend struct Foo <: AbstractFoo
    e::T
end

is really more like:

struct Foo{T} <: AbstractFoo
    a
    b::Int
    c::T
    d::Vector{T}
    e::T
end
=#
