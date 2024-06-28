#=
Created on Thr 14 Jun 2024
Updated on Thr 27 Jun 2024
Translated from Python (code dated 09/08/2017) and enhanced for Julia.
=#

"""
This module was taken from the author's course on numerical methods at TAMU.

From a biologic interpretation, the genetic algorithm implemented here is a
colony of creatures that advances their quality of life (fitness) from one
generation to the next.  Class 'GeneticAlgorithm,' defined at the end of this
file, is a colony or collection of creatures whose population is sustained from
one generation to the next.  Mating between creatures occurs through a process
known as tournament play, where the most fit contestant from a random selection
of contestants is chosen for mating.  Typically, each successive generation is
more fit than its predecessor, i.e., the colony's quality improves over time.

It is through class 'GeneticAlgorithm' that a user applies this genetic
algorithm, which is why this comment section has been listed at the top of this
file.  This colony of creatures has the following interface:

Constructor

    ga = GeneticAlgorithm(species, significantFigures, probabilityOfMutations,
                          probabilityOfCrossovers, probabilityOfImmigrants,
                          namesOfParameters, alienParameters, minParameters,
                          maxParameters)
    where
        species                 An implementation of the abstract class Species.
                                Here is where the user's model is introduced,
                                and where the experimental data are assigned.
        significantFigures      The number of significant figures of accuracy
                                sought in a final solution for the parameters.
                                There is a strong correlation between this value
                                and how long it takes to get a solution.  Values
                                of 4 or 5 are common.  They are bound to [1, 7].
        probabilityOfMutations  The probability of a gene mutation occurring at
                                conception, i.e. a gene swapping its expression.
                                Typically this is a low value, e.g., < 0.01.
        probabilityOfCrossovers The probability of a crossover event occurring
                                at conception, i.e., a chromosome splitting.
                                Typically this is large, e.g., > 0.85.
        probabilityOfImmigrants The probability of introducing an immigrant
                                into the gene pool after the first generation.
                                (The first generation is all immigrants.)  Use
                                as you see fit.  I usually assign it so as to
                                introduce one immigrant per generation or so.
        namesOfParameters       A string associated with each model parameter.
                                This makes the report much more readable.
        alienParameters         If the user has a best guess, then this is
                                where these parameters are input.  If the user
                                does not, then send: Vector{Float64}(undef,0).
        minParameters           Lower bounds for the parameters being sought.
        maxParameters           Upper bounds for the parameters being sought.

If a minParameters[i] equals its maxParameters[i], then parameter [i] is taken
to be fixed.

IMPORTANT!!!!!
For the parameter arrays, index 1 is reserved for selecting the best objective
function.  Specifically, a dimensionless error is computed via
    error[i] = (|model[i] - experiment[i]| / sdtDev(experiment))^p
where i indexes over the set of all experimental data being considered, with
parameter p (> 0) being a characteristic of the objective function.  It is
stored at index 1 in the various arrays for parameter data.  Once the error
vector is populated, then the mean error is calculated via
    meanError = mean(error)^(1/p)
thereby allowing for a quality value to be assigned as
    fitness = 1 / meanError
where the creature whose 'fitness' is the greatest is deemed the most fit.

Exponent 'p' is treated as another parameter to be optimized.  It can be fixed
or allowed to vary.  E.g., setting minParameters[1] = 2 and maxParameters[1] = 2
would impose a root mean squared error (RMSE) while setting minParameters[1] = 1
and maxParameters[1] = 1 would impose a linear error.  A more typical setting
would be minParameters[1] = 0.01 and maxParameters[1] = 100, which would seek
out an optimal objective function where the range of p lies within [0.01, 100]
for the specific set of data that are being considered.  I.e., it provides an
optimal objective function for a given model and set of data.

Methods

    p = bestCreatureParameters(ga)
        p       The full set of model parameters associated with the most fit
                or elite creature in the colony at its current generation.

    q = bestCreatureFitness(ga)
        q       The fitness parameter associated with the most fit or elite
                creature in the colony.  It is a measure of quality,
                specifically:
                    quality[i] = 1 / |modelResponse[i] - experimentResponse[i]|
                    quality    = sum over all quality[i]

    n = population(ga)
        n       The number of creatures that comprise a single generation.

    g = generation(ga)
        g       The number of generations that have been in existence.

    n = generationsToConvergence(ga)
        n       The expected number of generations needed for the genetic
                algorithm to converge upon a solution.

    advanceToNextGeneration!(ga)
                Advances the genetic algorithm to its next generation.  Giving
                the user access to this mechanism allows him/her a measure of
                control over how they choose to proceed, and what information
                they may wish to capture with each successive iteration in
                solution refinement.

    s = report(ga)
                This method returns a string 's' that comprises a report on the
                status of the optimizer at its current generation.  It should
                probably be called immediately after a colony object is created,
                and right after each advanceToNextGeneration call is made.  A
                wealth of statistics are reported here, specifically: optimal
                fitness, arithmetic mean, median, standard deviation, skewness
                and excess kurtosis, plus optimal parameters and their error
                bounds.  This provides the user with a means to make an
                assessment on the viability within a returned solution.
"""
module GeneticAlgorithms

using
    Statistics,
    StatsBase

import
    Printf: @sprintf

export
    # abstract type
    Species,

    # types
    Expression,
    Gene,
    Chromosome,
    Genome,
    Creature,
    GeneticAlgorithm,

    # operators
    ==, ≠,

    # external constructors
    procreate,
    alien,
    conceive,

    # methods
    copy,
    deepcopy,
    get,
    getindex,
    set!,
    setindex!,
    toString,

    isDominant,
    isRecessive,
    decode,
    encode!,
    mutate!,
    crossover,
    parameters,
    bestCreatureParameters,
    bestCreatrueFitness,
    population,
    generation,
    generationsToConvergence,
    advanceToNextGeneration!,
    report,

    # constants
    dominant,
    recessive

include("Expressions.jl")

include("Genes.jl")

include("Chromosomes.jl")

include("Genomes.jl")

include("Creatures.jl")

include("Species.jl")

# The genetic algorithm.

struct GeneticAlgorithm
    species::Species
    # fields pertaining to the colony
    significantFigures::Int64
    population::Int64
    generation::Int64
    # fields pertaining to the creatures of the colony
    elite::Creature
    adults::Vector{Creature}
    children::Vector{Creature}
    fitness::Vector{Float64}
    # fields pertaining to probabilities
    probabilityOfMutation::Float64
    probabilityOfCrossover::Float64
    probabilityOfImmigrant::Float64
    # fields pertaining to the parameters
    namesOfParameters::Vector{String}
    alienParameters::Vector{Float64}
    minParameters::Vector{Float64}
    maxParameters::Vector{Float64}

    # constructor

    function GeneticAlgorithm(species::Species, significantFigures::Int64, probabilityOfMutation::Float64, probabilityOfCrossover::Float64, probabilityOfImmigrant::Float64, namesOfParameters::Vector{String}, alienParameters::Vector{Float64}, minParameters::Vector{Float64}, maxParameters::Vector{Float64})

        # bound inputs

        if significantFigures < 1
            significantFigures = 1
        end
        if significantFigures > 7
            significantFigures = 7
        end

        if probabilityOfMutation < 1.0e-6
            probabilityOfMutation = 1.0e-6
        end
        if probabilityOfMutation > 0.999999
            probabilityOfMutation = 0.999999
        end

        if probabilityOfCrossover < 1.0e-6
            probabilityOfCrossover = 1.0e-6
        end
        if probabilityOfCrossover > 0.999999
            probabilityOfCrossover = 0.999999
        end

        if probabilityOfImmigrant < 1.0e-6
            probabilityOfImmigrant = 1.0e-6
        end
        if probabilityOfImmigrant > 0.999999
            probabilityOfImmigrant = 0.999999
        end

        # verify dimensions

        dimP = length(namesOfParameters)
        if (length(alienParameters) ≠ dimP) || (length(alienParameters) ≠ 0)
            msg = "The supplied alien has the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        if (length(minParameters) ≠ dimP) || (length(maxParameters) ≠ dimP)
            msg = "The supplied bounds have the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        # generation

        generation = 1

        # elite

        elite = alien(alienParameters, minParameters, maxParameters, significantFigures)

        # population

        alphabet = 2    # viz., dominant and recessive
        schemata = significantFigures
        # Formula is from D. Goldberg (2002) for estimating population size.
        population = Int64(ceil(alphabet^schemata * schemata * log(alphabet) + log(elite.genetics.genes)))

        # adults

        adults = Vector{Creature}(undef, population)
        adults[1] = elite
        for i in 2:population
            adult = procreate(minParameters, maxParameters, significantFigures)
            adult.fitness = _evaluate(species, adult)
            adults[i] = adult
            if adult.fitness > elite.fitness
                elite = adult
            end
        end

        # children

        children = Vector{Creature}(undef, population)

        # fitness

        fitness = Vector{Float64}(undef, population)
        for i in 1:population
            adult = adults[i]
            fitness[i] = adult.fitness
        end

        new(species, significantFigures, population, generation, elite, adults, children, fitness, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, namesOfParameters, alienParameters, minParameters, maxParameters)
    end
end # GeneticAlgorithm

# internal methods

function _tournamentPlay(ga::GeneticAlgorithm)::Creature
    combatants  = max(3, ga.population÷50)
    contestants = Vector{Creature}(undef, combatants)
    if ga.probabilityOfImmigrant > rand()
        mostFit = Creature(ga.minParameters, ga.maxParameters, ga.significantFigures)
        _evaluate(ga.species, mostFit)
    else
        for i in 1:combatants
            combatant = rand(1:ga.population)
            contestants[i] = ga.adults[combatant]
        end
        mostFit = contestants[1]
        for i in 2:combatants
            contestant = contestants[i]
            if contestant.fitness > mostFit.fitness
                mostFit = contestant
            end
        end
    end
    return mostFit
end # _tournamentPlay

function _mate!(ga::GeneticAlgorithm)
    # Elite creature from current generation lives into the next generation.
    ga.children[1] = deepcopy(ga.elite)
    ga.fitness[1]  = ga.elite.fitness
    for i in 2:population
        parentA = _tournamentPlay(ga)
        parentB = _tournamentPlay(ga)
        count = 0
        while parentB == parentA
            parentB = _tournamentPlay(ga)
            count = count + 1
            if count < (ga.population ÷ 2)
                println("The population is filling up with clones.")
                println("You should terminate and restart your run.")
                return nothing
            end
        end
        child = conceive(parentA, parentB, probabilityOfMutation, probabilityOfCrossover)
        child.fitness  = _evaluate(ga.species, child)
        ga.children[i] = child
        ga.fitness[i]  = child.fitness
        if child.fitness > ga.elite.fitness
            ga.elite = child
        end
    end
    return nothing
end # _mate!

function _evaluate(species::Species, c::Creature)::Float64
    par = parameters(c)
    runModel(species, par)
    if par[1] ≈ 0.0
        par[1] = sign(par[1]) * eps(Float32)
    end
    err = Vector{Matrix{Float64}}(undef, species.nExp)
    for i in 1:species.nExp
        err[i] = Matrix{Float64}(undef, species.nRsp[i], species.nPts[i])
    end
    n = 0
    quality = 0.0
    for i in 1:species.nExp
        for j in 1:species.nRsp[i]
            n = n + 1
            for k in 1:species.nPts[i]
                err[i][j,k] = ((abs(species.modR[i][j,k]
                    - species.expR[i][j,k]) / species.stdE[i][j])^par[1])
            end
            meanErr = mean(err[i][j,:])^(1.0/par[1])
            quality = quality + 1.0/meanErr
        end
    end
    return quality / n
end # _evaluate

function _stdDevElite(ga::GeneticAlgorithm)::Vector{Float64}
    bst = parameters(ga.elite)
    len = length(bst)
    sde = Vector{Float64}(undef, len)
    sum = zeros(Float64, len)
    for i in 1:ga.population
        adult = ga.adults[i]
        par_i = parameters(adult)
        for j in 1:len
            sum[j] = sum[j] + (par_i[j] - bst[j])^2
        end
    end
    for i in 1:len
        sde[i] = sqrt(sum[i] / (ga.population - 1))
    end
    return sde
end # _stdDevElite

function _2string(ga::GeneticAlgorithm, x::Float64)::String
    if (ga.significantFigures == 1) || (ga.significantFigures == 2)
        s = @sprintf "%.1e" x;
    elseif ga.significantFigures == 3
        s = @sprintf "%.2e" x;
    elseif ga.significantFigures == 4
        s = @sprintf "%.3e" x;
    elseif ga.significantFigures == 5
        s = @sprintf "%.4e" x;
    elseif ga.significantFigures == 6
        s = @sprintf "%.5e" x;
    else
        s = @sprintf "%.6e" x;
    end
    return s
end # _2string

# exported methods

function bestCreatureParameters(ga::GeneticAlgorithm)::Vector{Float64}
    return parameters(ga.elite)
end # bestCreatureParameters

function bestCreatrueFitness(ga::GeneticAlgorithm)::Float64
    return ga.elite.fitness
end # bestCreatureFitness

function population(ga::GeneticAlgorithm)::Int64
    return ga.population
end # population

function generation(ga::GeneticAlgorithm)::Int64
    return ga.generation
end # population

# Formula of D. Goldberg (2002) for estimating convergence.
function generationsToConvergence(ga::GeneticAlgorithm)::Int64
    combatants = max(3, ga.population÷50)
    convergeAt = Int64(ceil(sqrt(ga.elite.genetics.genes)
        * log(ga.population - combatants)))
    return convergeAt
end # generationsToConvergence

function advanceToNextGeneration!(ga::GeneticAlgorithm)
    ga.generation = ga.generation + 1
    _mate!(ga)
    for i in 1:population
        newAdult = ga.children[i]
        ga.adults[i] = copy(newAdult)
    end
end # advanceToNextGeneration

function report(ga::GeneticAlgorithm)::String
    s = "\n"
    s = string(s, "Statistics for generation ", ga.generation)
    s = string(s, " with a population size of ", ga.population, ".\n")
    s = string(s, "Optimal value and population statistics for fitness:\n")
    fitness = ga.elite.fitness
    if fitness ≥ 0.0
        s = string(s, "   optimal value    ", _2string(ga, fitness), "\n")
    else
        s = string(s, "   optimal value   ", _2string(ga, fitness), "\n")
    end
    statMean = mean(ga.fitness)
    if statMean ≥ 0.0
        s = string(s, "   arithmatic mean  ", _2string(ga, statMean), "\n")
    else
        s = string(s, "   arithmatic mean ", _2string(ga, statMean), "\n")
    end
    statMedian = median(ga.fitness)
    if statMedian ≥ 0.0
        s = string(s, "   median           ", _2string(ga, statMedian), "\n")
    else
        s = string(s, "   median          ", _2string(ga, statMedian), "\n")
    end
    statSkew = skewness(ga.fitness)
    if statSkew ≥ 0.0
        s = string(s, "   skewness         ", _2string(ga, statSkew), "\n")
    else
        s = string(s, "   skewness        ", _2string(ga, statSkew), "\n")
    end
    statKurtosis = kurtosis(ga.fitness)
    if statKurtosis ≥ 0.0
        s = string(s, "   excess kurtosis  ", _2string(ga, statKurtosis), "\n")
    else
        s = string(s, "   excess kurtosis ", _2string(ga, statKurtosis), "\n")
    end
    s = string(s, "Genetic genome from the most fit creature:\n")
    s = string(s, toString(ga.elite.genetics))
    s = string(s, "Listings of [minParameter, bestParameter, maxParameter]:\n")
    bestParameters = parameters(ga.elite)
    len = length(bestParameters)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", ga.namesOfParameters[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", ga.namesOfParameters[i], "\n")
            else
                s = string(s, i, ": ", ga.namesOfParameters[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, "[", _2string(ga.minParameters[i]), ", ",
            _2string(bestParameters[i]), ", ",
            _2string(ga.maxParameters[i]), "]\n")
    end
    s = string(s, "Values for bestParameter ± error, where an RMSE ")
    s = string(s, "is computed wrt best values.\n")
    s = string(s, "Data are the parameters from all living adults ")
    s = string(s, "in the current population.\n")
    err = _stdDevElite(ga)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", ga.namesOfParameters[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", ga.namesOfParameters[i], "\n")
            else
                s = string(s, i, ": ", ga.namesOfParameters[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, _2string(bestParameters[i]), " ± ",
            _2string(err[i]), "\n")
    end
    s = string(s, "\n")
    return s
end # report

end # module GeneticAlgorithms
