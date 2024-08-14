"""
struct ExperimentalData
    experiments::Int64                      # kinds/types of experimental data
    variables_control::Vector{Int64}        # control  variables per experiment
    variables_response::Vector{Int64}       # response variables per experiment
    data_points::Vector{Int64}              # experimental data  per experiment
    controls::Vector{Matrix{Float64}}       # controlled data:
                                                # [exp]x[nCtl[exp] x nPts[exp]]
    responses::Vector{Matrix{Float64}}      # response data:
                                                # [exp]x[nRes[exp] x nPts[exp]]
    responses_std::Vector{Vector{Float64}}  # standard deviation for responses
                                                # [exp]x[nRes[exp]]
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
    experiments::Int64                      # kinds of experimental data
    variables_control::Vector{Int64}        # control  variables per experiment
    variables_response::Vector{Int64}       # response variables per experiment
    data_points::Vector{Int64}              # experimental data  per experiment
    controls::Vector{Matrix{Float64}}       # controlled data:
                                                # [exp]x[nCtl[exp] x nPts[exp]]
    responses::Vector{Matrix{Float64}}      # response data:
                                                # [exp]x[nRes[exp] x nPts[exp]]
    responses_std::Vector{Vector{Float64}}  # standard deviation for responses
                                                # [exp]x[nRes[exp]]

    # constructor

    function ExperimentalData(experiments::Int64, variables_control::Vector{Int64}, variables_response::Vector{Int64}, data_points::Vector{Int64}, controls::Vector{Matrix{Float64}}, responses::Vector{Matrix{Float64}})

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

        responses_std = Vector{Vector{Float64}}(undef, experiments)
        for exp = 1:experiments
            responses_std[exp] = Vector{Float64}(undef, variables_response[exp])
            for res = 1:variables_response[exp]
                responses_std[exp][res] = std(responses[exp][res, :])
            end
        end

        new(experiments, variables_control, variables_response, data_points, controls, responses, responses_std)
    end

    function ExperimentalData(experiments::Integer, variables_control::Vector{Integer}, variables_response::Vector{Integer}, data_points::Vector{Integer}, controls::Vector{Matrix{Real}}, responses::Vector{Matrix{Real}})

        if (length(variables_control)     ≠ experiments
            || length(variables_response) ≠ experiments
            || length(data_points)        ≠ experiments
            || length(controls)           ≠ experiments
            || length(responses)          ≠ experiments)
            msg = "Vector lengths must equal the number of experiments."
            throw(DimensionMismatch, msg)
        end

        exp = convert(Int64, experiments)
        v_con = Vector{Int64}(undef, exp)
        v_res = Vector{Int64}(undef, exp)
        datPt = Vector{Int64}(undef, exp)
        for i in 1:exp
            v_con[i] = convert(Int64, variables_control[i])
            v_res[i] = convert(Int64, variables_response[i])
            datPt[i] = convert(Int64, data_points)
        end

        ctrl = Vector{Matrix}(undef, exp)
        for i in 1:exp
            matrix = Matrix{Float64}(undef, v_con[i], datPt[i])
            for j in 1:v_con[i]
                for k in 1:datPt[i]
                    matrix[j,k] = convert(Float64, controls[i][j,k])
                end
            end
            ctrl[i] = matrix
        end

        resp = Vector{Matrix}(undef, exp)
        for i in 1:exp
            matrix = Matrix{Float64}(undef, v_res[i], datPt[i])
            for j in 1:v_res[i]
                for k in 1:datPt[i]
                    matrix[j,k] = convert(Float64, responses[i][j,k])
                end
            end
            resp[i] = matrix
        end

        return ExperimentalData(exp, v_con, v_res, datPt, ctrl, resp)
    end

    function ExperimentalData(experiments::Int64, variables_control::Vector{Int64}, variables_response::Vector{Int64}, data_points::Vector{Int64}, controls::Vector{Matrix{Float64}}, responses::Vector{Matrix{Float64}}, responses_std::Vector{Vector{Float64}})

        new(exp, v_con, v_res, datPt, controls, responses, responses_std)::ExperimentalData
        end
    end
end # ExperimentalData

# Methods for storing and retrieving ExperimentalData to and from a file.

StructTypes.StructType(::Type{ExperimentalData}) = StructTypes.Struct()

"""
Method:\n
    toFile(data::GeneticAlgorithms.ExperimentalData, json_stream::IOStream)\n
writes a data structure `data` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    GeneticAlgorithms.toFile(data::GeneticAlgorithms.ExperimentalData, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
where <my_dir_path> is the path to your working directory wherein the file\n
<my_file_name> that is to be written to either exists or will be created,\n
and which must have a .json extension.
"""
function toFile(data::ExperimentalData, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, data)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

"""
Method:\n
    fromFile(::GeneticAlgorithms.ExperimentalData, json_stream::IOStream)\n
reads an instance of type Chromosome from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    data = GeneticAlgorithms.fromFile(::GeneticAlgorithms.ExperimentalData, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
which returns a `data,` an object of type GeneticAlgorithms.ExperimentalData.\n
Here <my_dir_path> is the path to your working directory wherein the file\n
to be read from, i.e., <my_file_name>, must exist, and which is to have a\n
.json extension.
"""
function fromFile(::Type{ExperimentalData}, json_stream::IOStream)::ExperimentalData
    if isopen(json_stream)
        data = JSON3.read(readline(json_stream), ExperimentalData)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return data
end

"""
AbstractParameters: a parent type for all Parameter types.
"""
abstract type AbstractParameters end

"""
Model is the type wherein an user establishes their model.
"""
struct Model{Parameters}
    θ::Parameters       # Must be mutable with Parameters <: AbstractParameters.
    d::ExperimentalData
end

#= 
It is vital that one's own Parameters type be mutable, e.g.,
    mutable struct MyParameters <: AbstractParameters
        a::Real
        b::Real
    end
from which a model is then created via
    mymodel = Model{MyParameters}(θ::MyParameters, d::ExperimentalData)
To assign or retrieve values as real vectors to/from a model's parameters, call
=#

function Base.:(get)(m::Model)::Vector{Float64}
    N = fieldcount(typeof(m.θ))
    θ = Vector{Float64}(undef, N)
    for n in 1:N
        symbol = fieldname(typeof(m.θ), n)
        θ[n]   = convert(Float64, getfield(m.θ, symbol))
    end
    return θ
end # get

function Base.:(getindex)(m::Model, index::Int)::Float64
    if index < 1 || index > fieldcount(typeof(m.θ))
        msg = "Index in getindex, i.e., in θ = m.θ[index], is out of range."
        throw(DimensionMismatch(msg))
    end
    symbol = fieldname(typeof(m.θ), index)
    θ = convert(Float64, getfield(m.θ, symbol))
    return θ
end # getindex

function set!(m::Model, θ::Vector{Float64})
    N = fieldcount(typeof(m.θ))
    if length(θ) == N
        for n in 1:N
            symbol = fieldname(typeof(m.θ), n)
            setfield!(m.θ, symbol, θ[n])
        end
    else
        msg = "The parameters to be assigned have the wrong dimension."
        error(msg)
    end
    return nothing
end # set!

function setindex!(m::Model, θ::Float64, index::Int)
    if index < 1 || index > fieldcount(typeof(m.θ))
        msg = "Index in setindex!, i.e., in m.θ[index] = θ, is out of range."
        throw(DimensionMismatch(msg))
    end
    symbol = fieldname(typeof(m.θ), index)
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
solve(m::Model) = solve(m.θ, m)

#=
As a template, consider:

function GeneticAlgorithms.solve(θ::MyParameters, m::Model)::Vector{Matrix{Float64}}
    # where MyParameters <: AbstractParameters
    # This function returns the model's response as a vector of matrices:
    #     [exp]x[nRes[exp] x nPts[exp]]  indexed as  [i][j,k]  where
    #         exp         experiment number with values  exp ∈ [1, nExp]
    #         nRes[exp]   number of  responses  for experiment exp
    #         nPts[exp]   number of data points for experiment exp

    # Create the array to be returned.

    responses_model = Vector{Matrix{Float64}}(undef, m.d.experiments)
    for exp in 1:m.d.experiments
        rows = m.d.variables_response[exp]
        cols = m.d.data_points[exp]
        response_matrix = Matrix{Float64}(undef, rows, cols)
        responses_model[exp] = response_matrix
    end

    # Solve your model using the parameters θ, subject to its m.d.controls.
    # This model can be a function, a differential equation, an integral
    # equation, or whatever else the user may choose to implement. Be aware that
    # here most likely lies the greatest expense of running a genetic algorithm.

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
The following was taken from url address: https://discourse.julialang.org/t/composition-and-inheritance-the-julian-way/11231/135

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
