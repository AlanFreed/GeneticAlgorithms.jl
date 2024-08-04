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
    3) conceived   (i.e., coming from its parents).
The first generation of a colony is procreated with the possible exception that
their is an alien (an Adam) assignment.  After that, all creatures are created
through conception, with a possibility of immigrants entering the population.

A creature has an interface of:

struct Creature
    fitness::Vector     # Measures of quality for the set of parameters, viz.,
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

    d = copy(c)         returns a copy 'd' of creature 'c'
    s = tostring(c)     returns a string 's' describing creature 'c'
    θ = phenotypes(c)   returns an array of the model's parameters 'θ'
"""
struct Creature
    fitness::Vector{Real}
    genetics::Genome

    # Internal constructors

    function Creature(parameters_min::Vector{Real}, parameters_max::Vector{Real}, significant_figures::Integer)

        fitness = Vector{Real}(undef, fitness_types)
        for i in 1:fitness_types
            fitness[i] = -1
        end

        genetics = Genome(parameters_min, parameters_max, significant_figures)

        new(fitness, genetics)
    end

    function Creature(fitness::Vector{Real}, genetics::Genome)
        new(fitness, genetics)
    end
end # Creature

# external constructors

function procreate(parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer)::Creature

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
            creature    = Creature(parameters_min, parameters_max, significant_figures)
            parameters_creature = phenotypes(creature)
            for i in 1:length(parameters_constrained)
                (pL, pR) = parameters_constrained[i]
                if parameters_creature[pL] > parameters_creature[pR]
                    recreate = true
                    break
                end
                if recreations == 25
                    # This is a safety valve.  It should not occur.
                    msg = "In 25 attempts, parameters θ[i] obeying "
                    msg = string(msg, "constraint θ[", pL, "]")
                    msg = string(msg, " < θ[", pR, "]\n were not gotten")
                    msg = string(msg, " through creature procreation.")
                    error(msg)
                end
            end
        end
    end

    return creature
end # procreate

function alien(parameters_alien::Vector{Real}, parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer)::Creature

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

function conceive(parentA::Creature, parentB::Creature, parameters_constrained::Vector{Tuple{Integer,Integer}}, probability_mutation::Real, probability_crossover::Real)::Creature

    genetics_child = crossover(parentA.genetics, parentB.genetics, probability_mutation, probability_crossover)

    fitness_child = Vector{Real}(undef, fitness_types)
    for i in 1:fitness_types
        fitness_child[i] = -1
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
            genetics_child = crossover(parentA.genetics, parentB.genetics, probability_mutation, probability_crossover)
            child = Creature(fitness_child, genetics_child)
            parameters_child = phenotypes(child)
            for i in 1:length(parameters_constrained)
                (pL, pR) = parameters_constrained[i]
                if parameters_child[pL] > parameters_child[pR]
                    reconceive = true
                    break
                end
                if reconceptions == 25
                    # This is a safety valve.  It should not occur.
                    msg = "In 25 attempts, parameters θ[i] obeying "
                    msg = string(msg, "constraint θ[", pL, "]")
                    msg = string(msg, " < θ[", pR, "]\n were not gotten")
                    msg = string(msg, " through creature conception.")
                    error(msg)
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
    fitness  = copy(c.fitness)
    genetics = copy(c.genetics)
    return Creature(fitness, genetics)
end # copy

function tostring(c::Creature)::String
    return tostring(c.genetics)
end # tostring

function phenotypes(c::Creature)::Vector{Real}
    θ = Vector{Real}(undef, c.genetics.chromosomes)
    for i in 1:c.genetics.chromosomes
        θ[i] = decode(c.genetics.genotypes[i])
    end
    return θ
end # phenotypes
