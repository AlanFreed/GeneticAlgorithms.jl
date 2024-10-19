"""
A creature has the genetic material of a genome, plus a means for its origin. After a creature (its object) has been created, it needs to have its genetic material assigned to it. This can occur in one of three ways:
1. procreation (i.e., random),
2. alien       (i.e., divinely assigned) or
3. conceived   (i.e., evolving from its parents).
The first generation of a colony is procreated with a possible exception that their may be an alien (an Adam) assignment. After that, all creatures are created through conception, with a possibility of an immigrant entering into the mating process.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# Creature

```julia
struct Creature
    fitness::PF.MReal
    DNA::Genome
end
```
where
1. *fitness* is the value of our objective function, which describes the goodness of fit as a sum of reciprocal norms, viz., 1/L₁ + 1/L₂ + 1/L∞, which are norms of the error between experimental response and model prediction.
2. *DNA* holds the genetic information of a creature.

## Constructors

### Internal constructors

```julia
creature = Creature(parameters_min::Vector{PF.PhysicalScalar},
                    parameters_max::Vector{PF.PhysicalScalar},
                    parameters_constrained::Vector{Tuple{Int,Int}}, 
                    significant_figures::Int)
```
```julia
creature = Creature(fitness::Float64, DNA::Genome)
```
```julia
creature = Creature(fitness::PF.MReal, DNA::Genome)
```
where
1. *parameters_min* is a vector containing the minimum scalar values for each phenotype (parameter) in the model being parameterized, i.e., θₘᵢₙ[i].
2. *parameters_max* is a vector containing the maximum scalar values for each phenotype (parameter) in the model being parameterized, i.e., θₘₐₓ[i].
3. *parameters_constrained* is a vector of tuples [(l,r), …] imposing constraints of θ[l] < θ[r], i.e., they are tuples of vector index pairs. If there are no such constraints for your model, send `parameters_constrained = Vector{Tuple}(undef, 0)`.
4. *significant_figures* specifies the number of digits in accuracy sought by a solution for the model's parameters.

> **NOTE**: if a `parameters_min[i] ≈ parameters_max[i]` then the phenotype being represented by the chromosome at index i will be taken to be constant, fixed to the value of this collapsed range.

### External constructors

```julia
creature = procreate(parameters_min::Vector{PF.PhysicalScalar},
                     parameters_max::Vector{PF.PhysicalScalar}, 
                     parameters_constrained::Vector{Tuple{Int,Int}}, 
                     significant_figures::Int)
```
> returns a creature of type Creature whose parameters are assigned at random. Its arguments are discussed above.

```julia
creature = alien(parameters_alien::Vector{PF.PhysicalScalar},
                 parameters_min::Vector{PF.PhysicalScalar},
                 parameters_max::Vector{PF.PhysicalScalar},
                 parameters_constrained::Vector{Tuple{Int,Int}},
                 significant_figures::Int)
```
> returns a creature of type Creature whose parameters are assigned by the user via vector `parameters_alien`. Its other arguments are discussed above.

```julia
creature = conceive(parentA::Creature, 
                    parentB::Creature,
                    parameters_constrained::Vector{Tuple{Int,Int}},
                    probability_mutation::Float64,
                    probability_crossover::Float64)
```
> returns a creature of type Creature whose parameters arise from a mating between two parents, viz., parentA and parentB. Here `probability_mutation` specifies the probability of a gene mutating at conception. This is typically quite small. While `probability_crossover` specifies the probability of the parents' chromosomes splitting and pairing at conception; hence, the child will share genes from both of its parents. This is typically quite large. Argument `parameters_constrained` is discussed above.

## Operators

`==` and `≠`

## Methods

```julia
cc = copy(creature::Creature)
```
> returns a copy cc of the creature.

```julia
str = toBinaryString(creature::Creature)
```
> returns a string str describing the creature written in a binary format.

```julia
str = toString(c::Creature)
```
> returns a string representation str of the phenotypes held by a creature.

```julia
θ = parameters(creature::Creature)
```
> returns an array of type Vector{PF.PhysicalScalar} for model parameters θ held in the DNA of a creature.

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

To write or read an instance of type *Creature* to or from a JSON file, call
```julia
toFile(creature, json_stream)
```
> which writes a creature of type *Creature* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
creature = fromFile(Creature, json_stream)
```
> reads a creature of type *Creature* from the JSON file attached to `json_stream`.
"""
struct Creature
    fitness::PF.MReal   # A measure of goodness of fit.
    DNA::Genome         # Complete set of genetic material.

    # Internal constructors

    function Creature(parameters_min::Vector{PF.PhysicalScalar},  
                      parameters_max::Vector{PF.PhysicalScalar}, 
                      parameters_constrained::Vector{Tuple{Int,Int}},
                      significant_figures::Int)
        fitness = PF.MReal()
        DNA = Genome(parameters_min, parameters_max,
                     parameters_constrained, significant_figures)
        new(fitness, DNA)::Creature
    end

    function Creature(fitness::Float64, DNA::Genome)
        new(PF.MReal(fitness), DNA)::Creature
    end
    
    function Creature(fitness::MReal, DNA::Genome)
        new(fitness, DNA)::Creature
    end
end # Creature

# external constructors

function procreate(parameters_min::Vector{PF.PhysicalScalar},     
                   parameters_max::Vector{PF.PhysicalScalar}, 
                   parameters_constrained::Vector{Tuple{Int,Int}}, 
                   significant_figures::Int)::Creature

    creature = Creature(parameters_min, parameters_max,
                        parameters_constrained, significant_figures)
    constraints = length(parameters_constrained)
    if constraints > 0
        pL = 0
        pR = 0
        recreate = false
        parameters_creature = parameters(creature)
        for i in 1:constraints
            (pL, pR) = parameters_constrained[i]
            if parameters_creature[pL] > parameters_creature[pR]
                recreate = true
                break
            end
        end
        recreations = 0
        while recreate
            recreate    = false
            recreations = recreations + 1
            if recreations == 25
                # This is a safety valve.  It should not occur in practice.
                msg = "In 25 attempts, parameters θ[i], which must obey "
                msg = string(msg, "constraint θ[", pL, "]")
                msg = string(msg, " < θ[", pR, "],\n were unable to be gotten")
                msg = string(msg, " through creature procreation.")
                error(msg)
            end
            creature = Creature(parameters_min, parameters_max,
                                parameters_constrained, significant_figures)
            parameters_creature = parameters(creature)
            for i in 1:constraints
                (pL, pR) = parameters_constrained[i]
                if parameters_creature[pL] > parameters_creature[pR]
                    recreate = true
                    break
                end
            end
        end
    end

    return creature
end # procreate

function alien(parameters_alien::Vector{PF.PhysicalScalar},
               parameters_min::Vector{PF.PhysicalScalar},
               parameters_max::Vector{PF.PhysicalScalar}, 
               parameters_constrained::Vector{Tuple{Int,Int}}, 
               significant_figures::Int)::Creature

    if length(parameters_alien) == 0
        creature = procreate(parameters_min, parameters_max,
                             parameters_constrained, significant_figures)
    else
        # Verify the input.
        if length(parameters_alien) ≠ length(parameters_min)
            msg = "Vector parameters_alien has the wrong length."
            throw(DimensionMismatch(msg))
        end

        # Ensure that the alien's parameters obey their constraints.
        for i in 1:length(parameters_constrained)
            (pL, pR) = parameters_constrained[i]
            if parameters_alien[pL] > parameters_alien[pR]
                msg = "Alien parameters θ[i] violate their constraint"
                msg = string(msg, " of θ[", pL, "] < θ[", pR, "].")
                error(msg)
            end
        end

        creature = Creature(parameters_min, parameters_max,
                            parameters_constrained, significant_figures)
        encode!(creature.DNA, parameters_alien)
    end

    return creature
end # alien

function conceive(parentA::Creature, 
                  parentB::Creature,
                  parameters_constrained::Vector{Tuple{Int,Int}},
                  probability_mutation::Float64,
                  probability_crossover::Float64)::Creature

    fitness_child = PF.MReal()
    DNA_child = crossover(parentA.DNA, parentB.DNA,
                          probability_mutation, probability_crossover)
    child = Creature(fitness_child, DNA_child)

    if length(parameters_constrained) > 0
        pL = 0
        pR = 0
        reconceive = false
        parameters_child = parameters(child)
        for i in 1:length(parameters_constrained)
            (pL, pR) = parameters_constrained[i]
            if parameters_child[pL] > parameters_child[pR]
                reconceive = true
                break
            end
         end
        reconceptions = 0
        while reconceive
            reconceive     = false
            reconceptions  = reconceptions + 1
            if reconceptions == 25
                # This is a safety valve.  It should not occur in practice.
                msg = "In 25 attempts, parameters θ[i], which must obey "
                msg = string(msg, "constraint θ[", pL, "]")
                msg = string(msg, " < θ[", pR, "],\n were not able to be ")
                msg = string(msg, " gotten through creature conception.")
                error(msg)
            end
            DNA_child = crossover(parentA.DNA, parentB.DNA,
                                  probability_mutation,
                                  probability_crossover)
            child = Creature(fitness_child, DNA_child)
            parameters_child = parameters(child)
            for i in 1:length(parameters_constrained)
                (pL, pR) = parameters_constrained[i]
                if parameters_child[pL] > parameters_child[pR]
                    reconceive = true
                    break
                end
            end
        end
    end

    return child
end # conceive

# operators

function Base.:(==)(creature_left::Creature, creature_right::Creature)::Bool
    if creature_left.DNA ≠ creature_right.DNA
        return false
    else
        return true
    end
end # ==

function Base.:≠(creature_left::Creature, creature_right::Creature)::Bool
    return !(creature_left == creature_right)
end # ≠

# methods

function Base.:(copy)(creature::Creature)::Creature
    fitness = copy(creature.fitness)
    DNA = copy(creature.DNA)
    return Creature(fitness, DNA)
end # copy

function toBinaryString(creature::Creature)::String
    return toBinaryString(creature.DNA)
end # toBinaryString

function toString(creature::Creature)::String
    return toString(creature.DNA)
end # toString

function parameters(creature::Creature)::Vector{PF.PhysicalScalar}
    return decode(creature.DNA)
end # parameters

# Methods for storing and retrieving a Creature to and from a file.

StructTypes.StructType(::Type{Creature}) = StructTypes.Struct()

function toFile(creature::Creature, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, creature)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{Creature}, json_stream::IOStream)::Creature
    if isopen(json_stream)
        creature = JSON3.read(readline(json_stream), Creature)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return creature
end

