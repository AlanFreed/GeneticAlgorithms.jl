"""
struct ExperimentalData
    experiments::Integer                    # kinds/types of experimental data
    variables_control::Vector{Integer}      # control  variables per experiment
    variables_response::Vector{Integer}     # response variables per experiment
    data_points::Vector{Integer}            # experimental data  per experiment
    controls::Vector{Matrix{Real}}          # controlled data:
                                            # [exp] x [nCtl[exp] x nPts[exp]]
    responses::Vector{Matrix{Real}}         # response data:
                                            # [exp] x [nRes[exp] x nPts[exp]]
    responses_std::Vector{Vector{Real}}     # standard deviation for responses
                                            # [exp] x [nRes[exp]]
end

Arrays variables_control, variables_response and data_points index as [i], where
    i ∈ [1, experiments].
Arrays controls and responses index as [i][j,k], where i ∈ [1, experiments],
    j ∈ [1, variables_control[i]] or [1, variables_response[i]], and where
    k ∈ [1, data_points[i]].
Array responses_std indexes as [i][j], where i ∈ [1, experiments] and where
    j ∈ [1, variables_response[i]].
"""
struct ExperimentalData
    experiments::Integer                    # kinds of experimental data
    variables_control::Vector{Integer}      # control  variables per experiment
    variables_response::Vector{Integer}     # response variables per experiment
    data_points::Vector{Integer}            # experimental data  per experiment
    controls::Vector{Matrix{Real}}          # controlled data:
                                            # [exp] x [nCtl[exp] x nPts[exp]]
    responses::Vector{Matrix{Real}}         # response data:
                                            # [exp] x [nRes[exp] x nPts[exp]]
    responses_std::Vector{Vector{Real}}     # standard deviation for responses
                                            # [exp] x [nRes[exp]]

    # constructor

    function ExperimentalData(experiments::Integer, variables_control::Vector{Integer}, variables_response::Vector{Integer}, data_points::Vector{Integer}, controls::Vector{Matrix{Real}}, responses::Vector{Matrix{Real}})

        # verify inputs

        if (length(variables_control)     ≠ experiments
            || length(variables_response) ≠ experiments
            || length(data_points)        ≠ experiments
            || length(controls)           ≠ experiments
            || length(responses)          ≠ experiments)
            msg = "Vector lengths must equal the number of experiments."
            throw(DimensionMismatch, msg)
        end
        for i in 1:experiments
            if size(controls[i]) ≠ (variables_control[i], data_points[i])
                msg = "Matrix size for controls is not compatible."
                throw(DimensionMismatch, msg)
            end
            if size(responses[i]) ≠ (variables_response[i], data_points[i])
                msg = "Matrix size for responses is not compatible."
                throw(DimensionMismatch, msg)
            end
        end

        # Create and populate the vector for standard deviations in response.

        responses_std = Vector{Vector{Real}}(undef, experiments)
        for exp = 1:experiments
            stdR[exp] = Vector{Real}(undef, variables_response[exp])
            for res = 1:variables_response[exp]
                responses_std[exp][res] = std(responses[exp][res, :])
            end
        end

        new(experiments, variables_control, variables_response, data_points, controls, responses, responses_std)
    end
end # ExperimentalData

"""
AbstractParameters: a parent type for all Parameter types.
"""
abstract type AbstractParameters end

"""
Model is the type wherein an user establishes their model.
"""
struct Model
    θ::Parameters <: AbstractParameters # This must be a mutable composite type.
    d::ExperimentalData
end

#= 
It is vital that one's own Parameters type be mutable, e.g.,
mutable struct MyParameters
    a::Real
    b::Real
end
along with an associated method
function Base.:(copy)(θ::MyParameters)::MyParameters
    ...
end # copy
from which a model is then created via
mymodel = Model(θ::MyParameters, d::ExperimentalData)
=#

function Base.:(get)(m::Model)::Vector{Real}
    N = fieldcount(m.θ)
    θ = Vector{Real}(undef, N)
    for n in 1:N
        symbol = fieldnames(m.θ, n)
        θ[n]   = getfield(m.θ, symbol)
    end
    return θ
end # get

function Base.:(getindex)(m::Model, index::Integer)::Real
    if index < 1 || index > fieldcount(m.θ)
        msg = "The index in getindex is out of range."
        throw(ErrorException, msg)
    end
    symbol = fieldnames(m.θ, index)
    θ = getfield(m.θ, symbol)
    return θ
end # getindex

function set!(m::Model, θ::Vector{Real})
    N = fieldcount(m.θ)
    if length(θ) == N
        for n in 1:N
            symbol = fieldnames(m.θ, n)
            setfield!(m.θ, symbol, θ[n])
        end
    else
        msg = "The parameters to be assigned have the wrong dimension."
        throw(ErrorException, msg)
    end
    return nothing
end # set!

function setindex!(m::Model, θ::Real, index::Integer)
    if index < 1 || index > fieldcount(m.θ)
        msg = "The index in setindex! is out of range."
        throw(ErrorException, msg)
    end
    symbol = fieldnames(m.θ, index)
    setfield!(m.θ, symbol, θ)
    return nothing
end # setindex!

"""
An user must create their own procedure that solves their model according to
the interface listed below.  See the examples directory for illustrations as to
how to accomplish this.  Multiple dispatch according to the arguments of 'solve'
determines which model is to be run.

The returned array for model response has the same dimensions as the array
m.d.responses, which holds the experimental response to be compared against.
"""
solve(m::Model) = solve(m.θ, m.d, m)::Vector{Matrix{Real}}

#=
As a template, consider:

function solve(m.θ::P, m.d::ExperimentalData, m::Model)::Vector{Matrix{Real}} where P <: AbstractParameters
    # This function returns the model's response as a vector of matrices:
    #   [exp] x [nRes[exp] x nPts[exp]]  indexed as  [i][j,k]  where
    #       exp         experiment number with values  exp ∈ [1, nExp]
    #       nRes[exp]   number of  responses  for experiment exp
    #       nPts[exp]   number of data points for experiment exp

    # Create this array.

    responses_model = Vector{Matrix{Real}}(undef, m.d.experiments)
    for exp in 1:m.d.experiments
        rows = m.d.variables_response[exp]
        cols = m.d.data_points[exp]
        response_matrix = Matrix{Real}(undef, rows, cols)
        responses_model[exp] = response_matrix
    end

    # Solve your model using the parameters θ, subject to its m.d.controls.
    # This model can be a function, a differential equation, an integral
    # equation, or whatever else the user may choose to implement.  Be aware
    # that here lies the greatest expense of running a genetic algorithm.

    # ...
    for i in 1:m.d.experiments
        # ...
        for j in 1:m.d.variables_response[i]
            # ...
            for k in 1:m.d.data_points[i]
                responses_model[i][j,k] = ...
            end
            # ...
        end
        # ...
    end
    # ...

    return responses_model
end # solve
=#

#=
The following was taken from the url address: https://discourse.julialang.org/t/composition-and-inheritance-the-julian-way/11231/135

Here's an example of how multiple dispatch interacts with inheritance of fields
and behaviors using composition.  Object has some fields and behaviors of its
own, and behaviors and even fields from Style objects.  These influence the
results of some_function using multiple dispatch.

abstract type Style end
struct StyleA <: Style end
struct StyleB <: Style
    val::Int
end

struct Object{S}
    style::S
    val::Bool
end

some_function(o::Object) = some_function(o.style, o)
some_function(::StyleA, o::Object) = o.val
some_function(t::StyleB, o) = t.val > 0

This function could be defined as a fallback for a new StyleC without any
changes to Object

some_function(c::StyleC, o) = c.val^2 > 10

And by using StyleC() as your style in Object will give you StyleC behaviors.
But you could additionally write a method to dispatch on both and handle the
interaction.

(if you're still wondering why, now think about adding another composed object
and dispatch for combinations of it and Style, and how you would just add
another argument to some_function and define a few methods, and it will work
and extend to whatever complexity of inheritance you need…)
=#
