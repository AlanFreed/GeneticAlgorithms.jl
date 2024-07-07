"""
This genetic algorithm calls internally the function:

function _evaluate(s::AbstractSpecies, c::Creature)
    ...
    θ = parameters(c)
    modelResponse = runModel(s, θ)
    ...
end

To use this genetic algorithm, a user must create their own data structure:

struct MySpecies <: AbstractSpecies
    nExp::Int64                     number of experiments
    nCtl::Vector{Int64}             number of control  variables per experiment
    nRes::Vector{Int64}             number of response variables per experiment
    nPts::Vector{Int64}             number of datum points per experiment
    ctrl::Vector{Matrix{Float64}}   control  data for  each experiment
    resp::Vector{Matrix{Float64}}   response data from each experiment
    stdR::Vector{Vector{Float64}}   standard deviation in the responses per exp
end

whose fields will be accessed by the genetic algorithm.

Arrays nCtl, nRes and nPts index as [i], where i ∈ [1, nExp].
Arrays ctrl and resp index as [i][j,k],  where i ∈ [1, nExp],
    j ∈ [1, nCtl nRes] and k ∈ [1, nPts].
Array stdR indexes as [i][j], where i ∈ [1, nExp] and j ∈ [1, nRes].

Along with the user's data structure, the user must supply their own function
(template as follows):

function runModel(s::MySpecies, θ::Vector{Float64})::Vector{Matrix{Float64}}
    # NOTE: the model's parameters are to index from 2,
    # as θ[1] is an internal parameter assigned by the optimizer.
    # Retrieve the data held within MySpecies s.
    nExp = experiments(s)
    nCtl = controlsPerExp(s)
    nRes = responsesPerExp(s)
    nPts = dataPointsPerExp(s)
    ctrl = controls(s)
    resp = responses(s)
    stdR = stdDevInResponses(s)
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
    # Run your model subject to MySpecies s' controls, viz., subject to ctrl.
    # This model can be a function, a differential equation, an integral
    # equation, or whatever else the user may choose to implement here.  Be
    # aware that here is where the greatest expense of this algorithm resides.
    # ...
    for i in 1:nExp
        # ...
        for j in 1:nRes[i]
            # ...
            for k in 1:nPts[i]
                # ...
            end
            # ...
        end
        # ...
    end
    # ...
    return modR
end
"""
abstract type AbstractSpecies end

# Dimension the problem to be solved.

experiments(s::AbstractSpecies) = s.nExp

controlsPerExp(s::AbstractSpecies) = s.nCtl

responsesPerExp(s::AbstractSpecies) = s.nRes

dataPointsPerExp(s::AbstractSpecies) = s.nPts

# Assign the experimental data.

controls(s::AbstractSpecies) = s.ctrl

responses(s::AbstractSpecies) = s.resp

stdDevInResponses(s::AbstractSpecies) = s.stdR

#= 
# The following template, of sorts, addresses how to create an optimizer using Julia's notion of multiple dispatch, versus an object oriented approach where the model is handled as a method.  This example was taken from the web site https://discourse.julialang.org/t/alternative-to-function-as-field-in-struct/55094/4

# If you have various models with different internal structures and different functions that can be applied to them, then you may still benefit from defining models and functions separately. You can use something like this

abstract type AbstractModel end

data(m::AbstractModel) = m.θ

# in optimize function
function optimize(m::AbstractModel)
    ...
    θ = data(m)
    # do something with θ
    ...
    # use function defined on `m`
    x = f(m)
end

# Of course, it would require a more involved setup for the concrete models. For each model, you have to define your own structure and necessary functions

struct MyModel{T} <: AbstractModel
    θ::Vector{T}
end

f(m::MyModel) = m.θ .^ 2

# But the benefit is that you will not be restricted by the initial type definition, i.e. you can add as many parameters as needed, and also there will be no issues with functions as a field.
=#
