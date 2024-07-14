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

    c = Creature(minParameters, maxParameters, constrainedParameters,
                 significantFigures, probabilityOfMutation,
                 probabilityOfCrossover)
or
    c = Creature(fitness, genetics)

External constructors

    c = procreate(minParameters, maxParameters, constrainedParameters,
                  significantFigures, probabilityOfMutation,
                  probabilityOfCrossover)
    c = alien(alienParameters, minParameters, maxParameters, 
              constrainedParameters, significantFigures, probabilityOfMutation,
              probabilityOfCrossover)
    c = conceive(parentA, parentB, constrainedParameters, probabilityOfMutation,
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
    fitness::Real
    genetics::Genome

    # Internal constructors

    function Creature(minParameters::Vector{Real}, maxParameters::Vector{Real}, significantFigures::Int)

        fitness  = -1.0
        genetics = Genome(minParameters, maxParameters, significantFigures)

        new(fitness, genetics)
    end

    function Creature(fitness::Real, genetics::Genome)

        new(fitness, genetics)
    end
end # Creature

# external constructors

function procreate(minParameters::Vector{Real}, maxParameters::Vector{Real}, constrainedParameters::Vector{Tuple{Int,Int}}, significantFigures::Int)::Creature

    creature = Creature(minParameters, maxParameters, significantFigures)

    if length(constrainedParameters) > 0
        recreate = false
        creatureParameters = parameters(creature)
        for i in 1:length(constrainedParameters)
            (pL, pR) = constrainedParameters[i]
            if creatureParameters[pL] > creatureParameters[pR]
                recreate = true
                break
            end
        end
        recreations = 1
        while recreate
            recreate    = false
            recreations = recreations + 1
            creature    = Creature(minParameters, maxParameters, significantFigures)
            creatureParameters = parameters(creature)
            for i in 1:length(constrainedParameters)
                (pL, pR) = constrainedParameters[i]
                if creatureParameters[pL] > creatureParameters[pR]
                    recreate = true
                    break
                end
                if recreations == 25
                    # This is a safety valve.  It should not occur.
                    msg = "In 25 attempts, parameters θ[i] obeying "
                    msg = string(msg, "constraint θ[", pL, "]")
                    msg = string(msg, " < θ[", pR, "]\n were not gotten")
                    msg = string(msg, " through creature procreation.")
                    throw(ExceptionError(msg))
                end
            end
        end
    end

    return creature
end # procreate

function alien(alienParameters::Vector{Real}, minParameters::Vector{Real}, maxParameters::Vector{Real}, constrainedParameters::Vector{Tuple{Int,Int}}, significantFigures::Int)::Creature

    if length(alienParameters) == 0
        creature = procreate(minParameters, maxParameters, constrainedParameters
            , significantFigures)
    else
        # Verify the input.
        if length(alienParameters) ≠ length(minParameters)
            msg = "Vector alienParameters has the wrong length."
            throw(DimensionMismatch, msg)
        end

        # Ensure that the alien's parameters obey their constraints.
        for i in 1:length(constrainedParameters)
            (pL, pR) = constrainedParameters[i]
            if alienParameters[pL] > alienParameters[pR]
                msg = "Alien parameters θ[i] violate their constraint"
                msg = string(msg, " of θ[", pL, "] < θ[", pR, "].")
                throw(ExceptionError(msg))
            end
        end

        creature = Creature(minParameters, maxParameters, significantFigures)
        encode!(creature.genetics, alienParameters)
    end

    return creature
end # alien

function conceive(parentA::Creature, parentB::Creature, constrainedParameters::Vector{Tuple{Int,Int}}, probabilityOfMutation::Real, probabilityOfCrossover::Real)::Creature

    fitness = -1.0

    childGenetics = crossover(parentA.genetics, parentB.genetics, probabilityOfMutation, probabilityOfCrossover)

    child = Creature(fitness, childGenetics)

    if length(constrainedParameters) > 0
        reconceive = false
        childParameters = parameters(child)
        for i in 1:length(constrainedParameters)
            (pL, pR) = constrainedParameters[i]
            if childParameters[pL] > childParameters[pR]
                reconceive = true
                break
            end
         end
        reconceptions = 1
        while reconceive
            reconceive    = false
            reconceptions = reconceptions + 1
            childGenetics = crossover(parentA.genetics, parentB.genetics, probabilityOfMutation, probabilityOfCrossover)
            child = Creature(fitness, childGenetics)
            childParameters = parameters(child)
            for i in 1:length(constrainedParameters)
                (pL, pR) = constrainedParameters[i]
                if childParameters[pL] > childParameters[pR]
                    reconceive = true
                    break
                end
                if reconceptions == 25
                    # This is a safety valve.  It should not occur.
                    msg = "In 25 attempts, parameters θ[i] obeying "
                    msg = string(msg, "constraint θ[", pL, "]")
                    msg = string(msg, " < θ[", pR, "]\n were not gotten")
                    msg = string(msg, " through creature conception.")
                    throw(ExceptionError(msg))
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

function parameters(c::Creature)::Vector{Real}
    phenotypes = Vector{Real}(undef, c.genetics.chromosomes)
    for i in 1:c.genetics.chromosomes
        phenotypes[i] = decode(c.genetics.genotypes[i])
    end
    return phenotypes
end # parameters
