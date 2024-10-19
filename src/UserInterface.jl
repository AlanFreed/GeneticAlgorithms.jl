"""
This data structure, for which the user is responsible, is where the actual data reside that are to be used to fit parameters belonging to the user's model.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# TheData

```julia
struct TheData
    experiments::Int
    conjugate_pairs::Vector{Int}
    data_points::Vector{Int}
    # independent variable
    time::Vector{PF.ArrayOfPhysicalScalars}
    # dependent variables
    control::Vector{Vector{PF.ArrayOfPhysicalScalars}}
    response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}}
    response_mod::Vector{Vector{PF.ArrayOfPhysicalScalars}}
end
```
where
1. *experiments* specifies the number of distinct experiments run. For example, one experiment could have been done in 1D, while a second experiment could have been done in 2D.
2. *conjugate_pairs* specifies the number of conjuage pairs, i.e., (input, output), (control, response), (cause, effect), … pairings that were monitored per experiment. Thermodynamic examples for conjugate pairs include (temperature, entropy), (pressure, volume), (chemical potential, particle number), etc.
3. *data_points* specifies the number of data points recorded per experiment over the conjugate pairs relevant to each experiment.
4. *time* is an ordered sequence of moments in time (the independent variable), one for each datum point in each experiment run.
5. *control* contains scalar fields (values with physical units) for the control variables of conjugate pairs at each datum point in each experiment run. These data will provide inputs to the user's model when it is being solved.
6. *response_exp* contains scalar fields (values with physical units) for the experimental response variables (subject to the control variables) of conjugate pairs at each datum point in each experiment run.
7. *response_mod* contains the model's response subject to the control variables. Their values are filled in by the solver whose histories, for many classes of model, are needed in order to advance a solution. These theoretical data will be contrasted with their experimental counterparts in order to assess a model's goodness of fit.

Arrays *conjugate_pairs* and *data_points* index as [i], where
> `i ∈ [1, experiments]`.

Array *time* indexes as [i][j], where
> `i ∈ [1, experiments]`, `j ∈ [1, data_points[i]]`.

Arrays *control*, *response_exp* and *response_mod* index as [i][j][k], where
> `i ∈ [1, experiments]`, `j ∈ [1, conjugate_pairs[i]]` and `k ∈ [1, data_points[i]]`, with `[i][j][:]` having the same physical units at each fixed i,j.

## Constructors

The user's constructor.
```julia
data = TheData(experiments::Int,
               conjugate_pairs::Vector{Int},
               data_points::Vector{Int},
               time::Vector{PF.ArrayOfPhysicalScalars},
               control::Vector{Vector{PF.ArrayOfPhysicalScalars}},
               response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}})
```

The constructor used, e.g., by JSON3.
```julia
data = TheData(experiments::Int,
               conjugate_pairs::Vector{Int},
               data_points::Vector{Int},
               time::Vector{PF.ArrayOfPhysicalScalars},
               control::Vector{Vector{PF.ArrayOfPhysicalScalars}},
               response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}},
               response_mod::Vector{Vector{PF.ArrayOfPhysicalScalars}})
```

## Operators

`==` and `≠`

## Methods

```julia
cc = copy(data::TheData)
```
> returns a copy cc of the data.

### Persistence

To open or close an IOStream attached to a JSON file, call
```julia
json_stream = PF.openJSONWriter(<my_dir_path>, <my_file_name.json>)
```
> which opens a `json_stream` of type *IOStream* for a file `<my_file_name.json>` located in directory `<my_dir_path>`, both of which are strings, while
```julia
PF.closeJSONStream(json_stream)
```
> flushes the buffer and closes this `json_stream`.

To write or read an instance of type *TheData* to or from a JSON file, call
```julia
toFile(data, json_stream)
```
> which writes the data of type *TheData* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
data = fromFile(TheData, json_stream)
```
> reads the data of type *TheData* from the JSON file attached to `json_stream`.
"""
struct TheData
    experiments::Int
    conjugate_pairs::Vector{Int}
    data_points::Vector{Int}
    # independent variable
    time::Vector{PF.ArrayOfPhysicalScalars}
    # dependent variables
    control::Vector{Vector{PF.ArrayOfPhysicalScalars}}
    response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}}
    response_mod::Vector{Vector{PF.ArrayOfPhysicalScalars}}

    # constructors

    function TheData(experiments::Int,
                     conjugate_pairs::Vector{Int},
                     data_points::Vector{Int},
                     time::Vector{PF.ArrayOfPhysicalScalars},
                     control::Vector{Vector{PF.ArrayOfPhysicalScalars}},
                     response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}})

        # Verify inputs.
        if experiments < 1
            error("Number of experiments must be positive.")
        end
        if (length(conjugate_pairs) ≠ experiments ||
            length(data_points)     ≠ experiments ||
            length(control)         ≠ experiments ||
            length(response_exp)    ≠ experiments)
            msg = "Vector lengths must equal the number of experiments."
            throw(DimensionMismatch(msg))
        end
        if length(time) ≠ 0 && length(time) ≠ experiments
            msg = "Length of vector time must be 0 (static case)\n or equal"
            msg = string(msg, " to the number of experiments (dynamic case).")
            throw(DimensionMismatch(msg))
        end
        for exp in 1:experiments
            if (length(control[exp])      ≠ conjugate_pairs[exp] ||
                length(response_exp[exp]) ≠ conjugate_pairs[exp])
                msg = "Number of independent and dependent variables "
                msg = string(msg, "were not equal for experiment ", exp, ".")
                throw(DimensionMismatch(msg))
            end
            for pair in 1:conjugate_pairs[exp]
                if length(time) == experiments
                    if (time[exp].array.len ≠ control[exp][pair].array.len ||
                        time[exp].array.len ≠ response_exp[exp][pair].array.len)
                        msg = "Time, control and response data were not paired "
                        msg = string(msg, "for experiment ", exp)
                        msg = string(msg, " at conjugate pair ", pair, ".")
                        throw(DimensionMismatch(msg))
                    end
                    for datum in 2:time[exp].array.len
                        if time[exp][datum] ≤ time[exp][datum-1]
                            error("Time must monotonically increase.")
                        end
                    end
                else
                    if (control[exp][pair].array.len ≠ 
                        response_exp[exp][pair].array.len)
                        msg = "Control and response data were not paired "
                        msg = string(msg, "for experiment ", exp)
                        msg = string(msg, " at conjugate pair ", pair, ".")
                        throw(DimensionMismatch(msg))
                    end
                end
            end
        end
        
        # Use the ξ statistic (from xicor.jl) to test for cause and effect.
        for exp in 1:experiments
            for pair in 1:conjugate_pairs[exp]
                x = control[exp][pair].array.vec
                y = response_exp[exp][pair].array.vec
                ξxy = Xicor.xicor(x, y)
                ξyx = Xicor.xicor(y, x)
                if ξyx > ξxy
                    println()
                    println("WARNING: According to Chatterjee's ξ statistic,")
                    println("your cause and effect are likely reversed for")
                    print("conjugate pair ", pair)
                    println(" in experiment ", exp, ".")
                end
            end
        end
        
        # Create a vector that will hold the model's responses.
        response_mod = Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef,
                              experiments)
        for exp = 1:experiments
            response_mod[exp] = Vector{PF.PhysicalScalar}(undef,
                                       conjugate_pairs[exp])
            pts = data_points[exp]
            for pair = 1:conjugate_pairs[exp]
                units = response_exp[exp][pair].units
                response_mod[exp][pair] = PF.ArrayOfPhysicalScalars(pts, units)
            end
        end

        new(experiments, conjugate_pairs, data_points,
            time, control, response_exp, response_mod)::TheData
    end

    function TheData(experiments::Int,
                     conjugate_pairs::Vector{Int},
                     data_points::Vector{Int},
                     time::Vector{PF.ArrayOfPhysicalScalars},
                     control::Vector{Vector{PF.ArrayOfPhysicalScalars}},
                     response_exp::Vector{Vector{PF.ArrayOfPhysicalScalars}},
                     response_mod::Vector{Vector{PF.ArrayOfPhysicalScalars}})

        new(experiments, conjugate_pairs, data_points,
            time, control, response_exp, response_mod)::TheData
    end
end # TheData

# operators

function Base.:(==)(data_left::TheData, data_right::TheData)::Bool
    if data_left.experiments ≠ data_right.experiments
        return false
    end
    for exp in 1:data_left.experiments
        if data_left.conjugate_pairs[exp] ≠ data_right.conjugate_pairs[exp]
            return false
        end
        if data_left.data_points[exp] ≠ data_right.data_points[exp]
            return false
        end
        if length(data_left.time) ≠ length(data_right.time)
            return false
        end
        if length(data_left.time) > 0
            for datum in 1:data_left.data_points[exp]
                if data_left.time[exp][datum] ≠ data_right.time[exp][datum]
                    return false
                end
            end
        end
        for pair in 1:data_left.conjugate_pairs[exp]
            for datum in 1:data_left.data_points[exp]
                if (data_left.control[exp][pair][datum] ≠
                    data_right.control[exp][pair][datum])
                    return false
                end
                if (data_left.response_exp[exp][pair][datum] ≠
                    data_right.response_exp[exp][pair][datum])
                    return false
                end
                if (data_left.response_mod[exp][pair][datum] ≠
                    data_right.response_mod[exp][pair][datum])
                    return false
                end
            end
        end
    end
    return true
end # ==

function Base.:≠(data_left::TheData, data_right::TheData)::Bool
    return !(data_left == data_right)
end # ≠

# methods
    
function Base.:(copy)(data::TheData)::TheData
    experiments     = copy(data.experiments)
    conjugate_pairs = Vector{Int}(undef, data.experiments)
    data_points     = Vector{Int}(undef, data.experiments)
    for exp in 1:data.experiments
        conjugate_pairs[exp] = copy(data.conjugate_pairs[exp])
        data_points[exp]     = copy(data.data_points[exp])
    end
    if length(data.time) == 0
        time = Vector{PF.ArrayOfPhysicalScalars}(undef, 0)
    else
        time = Vector{PF.ArrayOfPhysicalScalars}(undef, data.experiments)
    end
    control = (
        Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, data.experiments))
    response_exp = (
        Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, data.experiments))
    response_mod = (
        Vector{Vector{PF.ArrayOfPhysicalScalars}}(undef, data.experiments))
    for exp in 1:experiments
        if length(data.time) == data.experiments
            time[exp] = PF.ArrayOfPhysicalScalars(
                data.data_points[exp], data.time[exp].units)
        end
        control[exp] = (
            Vector{PF.ArrayOfPhysicalScalars}(undef, conjugate_pairs[exp]))
        response_exp[exp] = (
            Vector{PF.ArrayOfPhysicalScalars}(undef, conjugate_pairs[exp]))
        response_mod[exp] = (
            Vector{PF.ArrayOfPhysicalScalars}(undef, conjugate_pairs[exp]))
        for pair in 1:conjugate_pairs[exp]
            control[exp][pair] = PF.ArrayOfPhysicalScalars(
                data.data_points[exp], data.control[exp][pair].units)
            response_exp[exp][pair] = PF.ArrayOfPhysicalScalars(
                data.data_points[exp], data.response_exp[exp][pair].units)
            response_mod[exp][pair] = PF.ArrayOfPhysicalScalars(
                data.data_points[exp], data.response_mod[exp][pair].units)
            for datum in 1:data_points[exp]
                if length(data.time) == data.experiments
                    time[exp][datum] = PF.copy(data.time[exp][datum])
                end
                control[exp][pair][datum] = PF.copy(
                    data.control[exp][pair][datum])
                response_exp[exp][pair][datum] = PF.copy(
                    data.response_exp[exp][pair][datum])
                response_mod[exp][pair][datum] = PF.copy(
                    data.response_mod[exp][pair][datum])
            end
        end
    end
    thedata = TheData(experiments, conjugate_pairs, data_points,
                      time, control, response_exp, response_mod)
    return thedata
end # copy

# Methods for storing and retrieving TheData to and from a file.

StructTypes.StructType(::Type{TheData}) = StructTypes.Struct()

function toFile(data::TheData, json_stream::IOStream)
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

function fromFile(::Type{TheData}, json_stream::IOStream)::TheData
    if isopen(json_stream)
        data = JSON3.read(readline(json_stream), TheData)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return data
end

"""
A parent type for all parameter types that are to be supplied to the constructor of type Model.
"""
abstract type AbstractParameters end

"""
Here is where an user establishes their own model for parameterization.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF,
    GeneticAlgorithms as GA
```

# Model

```julia
struct Model{Parameters}
    θ::Parameters
    data::TheData
end
```
where
1. *θ* is an instance of type Parameters, which must be a subtype of supertype AbstractParameters. Each parameter in this struct is to be a scalar field of type PF.PhysicalScalar whose value is mutable but whose physical units are not.
2. *data* is an instance of type TheData, which is where the user supplies their data that is to be fit by this model::Model.

## Methods

```julia
θ = get(model::Model)
```
> returns a vector θ of type Vector{PF.PhysicalScalar} that holds the parameters belonging to a model, i.e., `θ = get(model)`.

```julia
θ = getindex(model::Model, index::Int)
```
> returns a scalar θ of type PF.PhysicalScalar held in the vector of parameters at location index in the model,  i.e., `θ = model[index]`.

```julia
set!(model::Model, θ::Vector{PF.PhysicalScalar})
```
> assigns a vector of parameters θ to a model, i.e., `model.θ = set!(model, θ)`.

```julia
setindex!(model::Model, θ::PF.PhysicalScalar, index::Int)
```
> assigns a scalar θ to location index in the vector of parameters belonging to a model, i.e., `model.θ[index] = θ`. 

### Solver

```julia
solve!(model::Model) = solve!(model.θ, model)
```
> The user must write their own function `solve!(model.θ, model)`. The compiler will use multiple dispatch on model.θ to select the correct implementation of solve at runtime.

> This method solves a model using the control data supplied within field model.data.control, placing the model's predicted response into field model.data.response_mod.

# Template

First, write a function 
```julia
function myData()::GA.TheData
``` 
that holds your experimental data. For example, `mydata = myData()` returned as an instance of type `GA.TheData`.

Second, create a type that holds your model's parameters, viz.,
```julia
struct MyParameters <: AbstractParameters
    # Enter parameters as instances of type `PF.PhysicalScalar`,
    # where the value of each scalar is mutable, but its units are not.
    …
end
```
Now, using your data (function `mydata = myData()`) and your model parameters (constructor `myparameters = MyParameters()`), create your model via
```julia
mymodel = GA.Model{MyParameters}(myparameters::MyParameters, mydata::GA.TheData)
```

At this point, one must write a solver for their model, which would be called as
```julia
GA.solve!(::MyParameters, mymodel::GA.Model)
```
> whose solution is placed into field `mymodel.data.response_mod`.

This solver should have a structure like:
```julia
function GA.solve!(::MyParameters, mymodel::GA.Model)
    # Type MyParameters allows Julia to select solve! via multiple dispatch.
    for exp in 1:mymodel.data.experiments
        for pair in 1:mymodel.data.conjugate_pairs[exp]
            for datum in 1:mymodel.data.data_points[exp]
                y = _mysolve(mymodel; exp, pair, datum)
                mymodel.data.response_mod[exp][pair][datum] = y
            end
        end
    end
    return nothing
end # GA.solve!
```
> where it is important that one creates a function `GA.solve!`, not just `solve!`.

This general solver calls a model specific solver, e.g.,
```julia
function _mysolver(mymodel::GA.Model; 
                   exp::Int, pair::Int, datum::Int)::PF.PhysicalScalar
    # Write your solver for 
    #    y = f(mymodel; exp, pair, datum)
    # where f is what you're writing here, which may solve a function, 
    # a differential equation, an integral equation, etc., given arguments:
    # 1. mymodel is the user's model,
    # 2. exp denotes the experiment,
    # 3. pair denotes the conjugate pair of this experiment, and
    # 4. datum specifies the data point along its solution path.
    # With this information, the user can writes a solver that solves for
    # 5. y is the response or output predicted by the model.
    # A model is said to be static if it does not depend upon time; 
    # otherwise, the model is said to be dynamic, with time being
    # stored in the field: mymodel.data.time.
    y = …
    return y
end # _mysolver
```
"""
struct Model{Parameters}
    θ::Parameters # Parameters <: AbstractParameters
    data::TheData
end

# Methods

function Base.:(get)(model::Model)::Vector{PF.PhysicalScalar}
    N = fieldcount(typeof(model.θ))
    θ = Vector{PF.PhysicalScalar}(undef, N)
    for n in 1:N
        symbol = fieldname(typeof(model.θ), n)
        θ[n]   = getfield(model.θ, symbol)
    end
    return θ
end # get

function Base.:(getindex)(model::Model, index::Int)::PF.PhysicalScalar
    if index < 1 || index > fieldcount(typeof(model.θ))
        msg = "Index in getindex, i.e., in θ = model.θ[index], is out of range."
        throw(DimensionMismatch(msg))
    end
    symbol = fieldname(typeof(model.θ), index)
    θ = getfield(model.θ, symbol)
    return θ
end # getindex

function set!(model::Model, θ::Vector{PF.PhysicalScalar})
    N = fieldcount(typeof(model.θ))
    if length(θ) == N
        for n in 1:N
            symbol = fieldname(typeof(model.θ), n)
            setfield!(model.θ, symbol, θ[n])
        end
    else
        error("The parameters to be assigned have the wrong dimension.")
    end
    return nothing
end # set!

function setindex!(model::Model, θ::PF.PhysicalScalar, index::Int)
    if index < 1 || index > fieldcount(typeof(model.θ))
        msg = "Index in setindex!, i.e., "
        msg = string(msg, "in model.θ[index] = θ, is out of range.")
        throw(DimensionMismatch(msg))
    end
    symbol = fieldname(typeof(model.θ), index)
    setfield!(model.θ, symbol, θ)
    return nothing
end # setindex!

# The user must create their own procedure that solves their model according
# to the interface listed below, i.e., the user must write their own function 
# solve!(model.θ, model). See the examples directory for illustrations as to
# how to accomplish this. Multiple dispatch, according to the arguments of
# method solve!, determines which model is to be run.

solve!(model::Model) = solve!(model.θ, model)

#=
The following commentary was taken from url address:
https://discourse.julialang.org/t/composition-and-inheritance-the-julian-way/11231/135
It is useful in helping to understand the logic that is being used here.

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
