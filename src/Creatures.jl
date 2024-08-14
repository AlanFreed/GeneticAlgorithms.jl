"""
There are four possible measures of fitness, i.e., objective functions, from
which one can choose, with the default being set at 3 when creating a colony.
    1) Minimize expectation in magnitude of error.
    2) Minimize expectation in squared error.
    3) Minimize variance in error.
    4) Maximize covariance between experimental response and model prediction.
Values assigned to the fitness vector are evaluated during the mating process.
"""
const fitness_types = 4

"""
A creature has the genetic material of a genome, plus a means for its origin.
After a creature (its object) has been created, it needs to have its genetic
material assigned to it.  This can occur in one of three ways:
    1) procreation (i.e., random),
    2) alien       (i.e., divinely assigned) or
    3) conceived   (i.e., evolving from its parents).
The first generation of a colony is procreated with the possible exception that
their is an alien (an Adam) assignment.  After that, all creatures are created
through conception, with a possibility of immigrants entering the population.

A creature has an interface of:

struct Creature
    fitness::MVector    # Measures of quality for the set of parameters, viz.,
                        # fitness = [ϕ₁ ϕ₂ ϕ₃ ϕ₄] where:
                        #    ϕ₁ = 1 / expectation for the magnitude of error,
                        #    ϕ₂ = 1 / expectation for the squared error,
                        #    ϕ₃ = 1 / variance in the error,
                        #    ϕ₄ = the covariance between experiment and model.
    genetics::Genome    # The genetic information: an array of chromosomes.
end

Internal constructor

    c = Creature(parameters_min, parameters_max, parameters_constrained,
                 significant_figures, probability_mutation,
                 probability_crossover)

where  (a, b) = parameters_constrained[i]  imposes the constraint  θ[a] < θ[b]

External constructors

    c = procreate(parameters_min, parameters_max, parameters_constrained,
                  significant_figures, probability_mutation,
                  probability_crossover)
    c = alien(parameters_alien, parameters_min, parameters_max, 
              parameters_constrained, significant_figures, probability_mutation,
              probability_crossover)
    c = conceive(parentA, parentB, parameters_constrained, probability_mutation,
                 probability_crossover)

Operators

    ==, ≠

Methods

    d = copy(c)             returns a copy 'd' of creature 'c'
    s = toBinaryString(c)   returns a string 's' describing creature 'c'
    θ = phenotypes(c)       returns an array of the model's parameters 'θ'
"""
struct Creature
    fitness::PhysicalFields.MVector     # A mutable vector of Float64 elements.
    genetics::Genome

    # Internal constructors

    function Creature(parameters_min::Vector{Float64}, parameters_max::Vector{Float64}, significant_figures::Int64)

        fitness = PhysicalFields.MVector(fitness_types)
        for i in 1:fitness_types
            fitness[i] = convert(Float64, -1)   # -1 implies not yet assigned
        end

        genetics = Genome(parameters_min, parameters_max, significant_figures)

        new(fitness, genetics)::Creature
    end

    function Creature(parameters_min::Vector{Real}, parameters_max::Vector{Real}, significant_figures::Integer)

        if length(parameters_min) ≠ length(parameters_max)
            msg = "Vectors parameters_min and parameters_max must have the same length."
            throw(DimensionMismatch, msg)
        end

        parameters = length(parameters_min)
        p_min = Vector{Float64}(undef, parameters)
        p_max = Vector{Float64}(undef, parameters)
        for i in 1:parameters
            p_min[i] = convert(Float64, parameters_min[i])
            p_max[i] = convert(Float64, parameters_max[i])
        end
        sigfig = convert(Int64, significant_figures)

        return Creature(p_min, p_max, sigfig)
    end

    function Creature(fitness::MVector, genetics::Genome)

        if length(fitness) ≠ fitness_types
            msg = "Argument fitness must have a length of "
            msg = string(msg, fitness_types, ".")
            error(msg)
        end

        new(fitness, genetics)::Creature
    end
end # Creature

# external constructors

function procreate(parameters_min::Vector{Float64}, parameters_max::Vector{Float64}, parameters_constrained::Vector{Tuple{Int64,Int64}}, significant_figures::Int64)::Creature

    creature = Creature(parameters_min, parameters_max, significant_figures)

    if length(parameters_constrained) > 0
        recreate = false
        parameters_creature = phenotypes(creature)
        for i in 1:length(parameters_constrained)
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
                # This is a safety valve.  It should not occur.
                msg = "In 25 attempts, parameters θ[i] obeying "
                msg = string(msg, "constraint θ[", pL, "]")
                msg = string(msg, " < θ[", pR, "]\n were not gotten")
                msg = string(msg, " through creature procreation.")
                error(msg)
            end
            creature = Creature(parameters_min, parameters_max, significant_figures)
            parameters_creature = phenotypes(creature)
            for i in 1:length(parameters_constrained)
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

function alien(parameters_alien::Vector{Float64}, parameters_min::Vector{Float64}, parameters_max::Vector{Float64}, parameters_constrained::Vector{Tuple{Int64,Int64}}, significant_figures::Int64)::Creature

    if length(parameters_alien) == 0
        creature = procreate(parameters_min, parameters_max, parameters_constrained, significant_figures)
    else
        # Verify the input.
        if length(parameters_alien) ≠ length(parameters_min)
            msg = "Vector parameters_alien has the wrong length."
            throw(DimensionMismatch, msg)
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

        creature = Creature(parameters_min, parameters_max, significant_figures)
        encode!(creature.genetics, parameters_alien)
    end

    return creature
end # alien

function conceive(parentA::Creature, parentB::Creature, parameters_constrained::Vector{Tuple{Int64,Int64}}, probability_mutation::Float64, probability_crossover::Float64)::Creature

    genetics_child = crossover(parentA.genetics, parentB.genetics, probability_mutation, probability_crossover)

    fitness_child = PhysicalFields.MVector(fitness_types)
    for i in 1:fitness_types
        fitness_child[i] = convert(Float64, -1) # -1 implies not yet assigned
    end

    child = Creature(fitness_child, genetics_child)

    if length(parameters_constrained) > 0
        reconceive = false
        parameters_child = phenotypes(child)
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
                # This is a safety valve.  It should not occur.
                msg = "In 25 attempts, parameters θ[i] obeying "
                msg = string(msg, "constraint θ[", pL, "]")
                msg = string(msg, " < θ[", pR, "]\n were not gotten")
                msg = string(msg, " through creature conception.")
                error(msg)
            end
            genetics_child = crossover(parentA.genetics, parentB.genetics, probability_mutation, probability_crossover)
            child = Creature(fitness_child, genetics_child)
            parameters_child = phenotypes(child)
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

function Base.:(==)(cL::Creature, cR::Creature)::Bool
    if cL.genetics ≠ cR.genetics
        return false
    else
        return true
    end
end # ==

function Base.:≠(cL::Creature, cR::Creature)::Bool
    if cL.genetics ≠ cR.genetics
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(copy)(c::Creature)::Creature
    fitness  = PhysicalFields.copy(c.fitness)
    genetics = copy(c.genetics)
    return Creature(fitness, genetics)
end # copy

function toBinaryString(c::Creature)::String
    return toBinaryString(c.genetics)
end # toBinaryString

function phenotypes(c::Creature)::Vector{Float64}
    θ = Vector{Float64}(undef, c.genetics.chromosomes)
    for i in 1:c.genetics.chromosomes
        θ[i] = decode(c.genetics.genotypes[i])
    end
    return θ
end # phenotypes

# Methods for storing and retrieving a Creature to and from a file.

StructTypes.StructType(::Type{Creature}) = StructTypes.Struct()

"""
Method:\n
    toFile(creature::GeneticAlgorithms.Creature, json_stream::IOStream)\n
writes a data structure `creature` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    GeneticAlgorithms.toFile(creature::GeneticAlgorithms.Creature, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
where <my_dir_path> is the path to your working directory wherein the file\n
<my_file_name> that is to be written to either exists or will be created,\n
and which must have a .json extension.
"""
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

"""
Method:\n
    fromFile(c::GeneticAlgorithms.Creature, json_stream::IOStream)\n
reads an instance of type Creature from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>::String, <my_file_name>::String)\n
    ...\n
    creature = GeneticAlgorithms.fromFile(::GeneticAlgorithms.Creature, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
which returns a `creature,` an object of type GeneticAlgorithms.Creature.\n
Here <my_dir_path> is the path to your working directory wherein the file\n
to be read from, i.e., <my_file_name>, must exist, and which is to have a\n
.json extension.
"""
function fromFile(::Type{Creature}, json_stream::IOStream)::Creature
    if isopen(json_stream)
        creature = JSON3.read(readline(json_stream), Creature)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return creature
end
