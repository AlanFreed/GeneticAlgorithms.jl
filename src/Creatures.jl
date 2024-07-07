"""
A creature has the genetic material of a genome plus a means for its origin.
After a creature (its object) has been created, it needs to have its genetic
material assigned to it.  This can occur in one of three ways:
    1) procreation (i.e., random),
    2) alien       (i.e., divinely assigned) or
    3) conceived   (i.e., coming from its parents).
The first generation of a colony is procreated with the possible exception that
their be an alien (an Adam) assignment.  After that, all creatures are created
through conception.

A creature has an interface of:

mutable struct Creature
    fitness             # a measure of quality for the set of parameters
    genetics::Genome    # the genetic information: an array of chromosomes
end

Internal constructors

    c = Creature(minParameters, maxParameters, significantFigures,
                 probabilityOfMutation, probabilityOfCrossover)
or
    c = Creature(fitness, genetics)

External constructors

    c = procreate(minParameters, maxParameters, significantFigures,
                  probabilityOfMutation, probabilityOfCrossover)
    c = alien(alienParameters, minParameters, maxParameters,
              significantFigures, probabilityOfMutation, probabilityOfCrossover)
    c = conceive(parentA, parentB, probabilityOfMutation,
                 probabilityOfCrossover)

Operators

    ==, ≠

Methods

    d = copy(c)         returns a copy 'd' of creature 'c'
    d = deepcopy(c)     returns a deep copy 'd' of creature 'c'
    s = toString(c)     returns a string 's' describing creature 'c'
    θ = parameters(c)   returns an array of all parameters 'θ'
"""
mutable struct Creature
    fitness::Float64
    genetics::Genome

    # Internal constructors

    function Creature(minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64)

        fitness  = -1.0
        genetics = Genome(minParameters, maxParameters, significantFigures)

        new(fitness, genetics)
    end

    function Creature(fitness::Float64, genetics::Genome)

        new(fitness, genetics)
    end
end # Creature

# external constructors

function procreate(minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64)::Creature

    creature = Creature(minParameters, maxParameters, significantFigures)

    return creature
end # procreate

function alien(alienParameters::Vector{Float64}, minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64)::Creature

    if length(alienParameters) ≠ length(minParameters)
        msg = "Vector alienParameters has the wrong length."
        throw(DimensionMismatch, msg)
    end

    creature = Creature(minParameters, maxParameters, significantFigures)

    if length(alienParameters) ≠ 0
        encode!(creature.genetics, alienParameters)
    end

    return creature
end # alien

function conceive(parentA::Creature, parentB::Creature, probabilityOfMutation::Float64, probabilityOfCrossover::Float64)::Creature

    fitness = -1.0

    childGenetics = crossover(parentA.genetics, parentB.genetics, probabilityOfMutation, probabilityOfCrossover)

    child = Creature(fitness, childGenetics)

    return child
end # conceive

# operators

function Base.:(==)(cL::Creature, cR::Creature)::Bool
    if cL.genetics ≠ cR.genetics
        return false
    end
    return true
end # ==

function Base.:≠(cL::Creature, cR::Creature)::Bool
    if !(cL == cR)
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(copy)(c::Creature)::Creature
    fitness  = copy(c.fitness)
    genetics = copy(c.genetics)
    return Creature(fitness, genetics)
end # copy

function Base.:(deepcopy)(c::Creature)::Creature
    fitness  = deepcopy(c.fitness)
    genetics = deeepcopy(c.genetics)
    return Creature(fitness, genetics)
end # deepcopy

function toString(c::Creature)::String
    return toString(c.genetics)
end # toString

function parameters(c::Creature)::Vector{Float64}
    phenotypes = Vector{Float64}(undef, c.genetics.chromosomes)
    for i in 1:c.genetics.chromosomes
        phenotypes[i] = decode(c.genetics.genotypes[i])
    end
    return phenotypes
end # parameters
