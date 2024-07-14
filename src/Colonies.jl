"""
A colony is a collection of creatures in a genetic algorithm.  It is the colony
that evolves from generation to generation, leading toward an optimum condition.

A Colony is a mutable data structure, because its data change from generation to
generation.  A colony has an interface of:

mutable struct Colony
    # Fields pertaining to the creatures of a colony.
    species                 # kind of creatures that comprise the colony
    probabilityOfMutation   # probability of single gene mutating at creation
    probabilityOfCrossover  # probability of a child sharing its parents genes
    probabilityOfImmigrant  # probability of an immigrant entering the colony
    parameterNames          # string representation for each model parameter
    minParameters           # lower bounds when searching for optimal parameters
    maxParameters           # upper bounds when searching for optimal parameters
    constrainedParameters   # tuples of specifying constraints: left < right
    significantFigures      # digits of accuracy sought for optimal parameters
    # Fields pertaining to the colony itself.
    populationSize          # number of creatures comprising the colony
    generation              # counter specifying the current generation
    generationsToConvergence    # generation whose elite is the solution sought
    elite                   # most fit resident in the current generation
    children                # adults for the next generation
    adults                  # adults from the current generation
    fitnesses               # fitness of each adult in the current generation
end

Constructors

    The constructor most likely to be called by a user.

    c = Colony(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, alienParameters, minParameters, maxParameters, constrainedParameters, significantFigures)

    where

        species                 An implementation of the class AbstractSpecies.
                                Here is where the user's model is introduced,
                                and where the experimental data are assigned.
        probabilityOfMutation   The probability of a gene mutation occurring at
                                conception, i.e. a gene swapping its expression.
                                Typically this is a low value, e.g., < 0.01.
        probabilityOfCrossover  The probability of a crossover event occurring
                                at conception, i.e., a chromosome splitting.
                                Typically this is a large value, e.g., > 0.85.
        probabilityOfImmigrant  The probability of introducing an immigrant
                                into the gene pool after the first generation.
                                (The first generation is all immigrants.)  Use
                                as you see fit.  I usually assign it so as to
                                introduce one immigrant per generation or so.
        parameterNames          A string associated with each model parameter.
                                This makes the report much more readable.
        alienParameters         If the user has a best guess for θ, then this is
                                where these parameters are input.  If the user
                                doesn't, then send: θ = Vector{Real}(undef,0).
        minParameters           Lower bounds for the parameters θ being sought.
        maxParameters           Upper bounds for the parameters θ being sought.
        constrainedParameters   Tuples of indices (left, right) that impose an
                                inequality constraint, viz.: θ[left] < θ[right].
        significantFigures      The number of significant figures of accuracy
                                sought in a final solution for the parameters.
                                There is a strong correlation between this value
                                and how long it takes to get a solution.  Values
                                of 4 or 5 are common.  They are bound to [1, 7].

    If parameter minParameters[i] equals parameter maxParameters[i] for any i, then parameter θ[i] is taken to be fixed.

The basic constructor (which will typically not be called externally) is

    c = Colony(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, constrainedParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitnesses)

IMPORTANT!!!!!

For the parameter arrays, index 1 is reserved for selecting the best objective
function.  Specifically, a dimensionless error is computed as
    error[i] = |model[i] - experiment[i]| / sdtDev(experiment)
where i indexes over the set of all experimental data being considered, whose
error is normalized by the standard deviation of experimental response.
A p-norm of this error is then computed, i.e.,
    ‖x‖ₚ = (∑ᵢ₌₁ⁿ|xᵢ|^p)^(1/p), p ≥ 1, with xᵢ = error[i] and n = s.nPts
from which a quality value is tabulated as
    fitness = ∑ᵢ₌₁ⁿ 1/‖x‖ₚ[i] with n = s.nExp
where the creature whose 'fitness' is the greatest is deemed most fit, i.e.,
the elite.  The first entry in all parameter arrays is a 'p' for the p-norm.

Exponent 'p', p ≥ 1, is treated here as just another parameter to be optimized.
It can be fixed or allowed to vary.  E.g., setting minParameters[1] = 1
and maxParameters[1] = 1 would impose the taxicab norm, while setting
minParameters[1] = 2 and maxParameters[1] = 2 would impose the Euclidean norm.
In addition, as 'p' goes toward infinity, the resulting p-norm will tend
towards the infinity (or maximum) norm.  A setting of minParameters[1] = 1 and
maxParameters[1] = 100 would be fairly typical.  Here one would seek an optimal
objective function (i.e., reciprocal p-norm) wherein the range of p would lie 
within the range [1, 100] for the specific set of data being considered.

Methods

advanceToNextGeneration!(c)     advances the colony 'c' to its next generation

str = report(c)                 returns a human-readable report via a string
                                'str,' which describes health of the current
                                generation and its statistics
"""
mutable struct Colony
    # Fields pertaining to the creatures of a colony.
    species::AbstractSpecies
    probabilityOfMutation::Real
    probabilityOfCrossover::Real
    probabilityOfImmigrant::Real
    parameterNames::Vector{String}
    minParameters::Vector{Real}
    maxParameters::Vector{Real}
    constrainedParameters::Vector{Tuple{Int,Int}}
    significantFigures::Int
    # Fields pertaining to the colony itself.
    populationSize::Int
    generation::Int
    generationsToConvergence::Int
    elite::Creature
    children::Vector{Creature}
    adults::Vector{Creature}
    fitnesses::Vector{Real}

    function Colony(species::AbstractSpecies, probabilityOfMutation::Real, probabilityOfCrossover::Real, probabilityOfImmigrant::Real, parameterNames::Vector{String}, alienParameters::Vector{Real}, minParameters::Vector{Real}, maxParameters::Vector{Real}, constrainedParameters::Vector{Tuple{Int,Int}}, significantFigures::Int)

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

        # Verify dimensions.
        # Vector lengths should be 1 + number of parameters in the model itself.

        dim = length(parameterNames)
        str = "p in p-norm, i.e., ‖x‖ₚ = (∑ᵢ₌₁ⁿ|xᵢ|^p)^(1/p),  p ≥ 1"
        parameterNames[1] = str

        if (length(alienParameters) ≠ dim) && (length(alienParameters) ≠ 0)
            msg = "The supplied alien has the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        if (length(minParameters) ≠ dim) || (length(maxParameters) ≠ dim)
            msg = "The supplied bounds have the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        # Ensure the p-norm has a value of p ≥ 1.
        if minParameters[1] < 1.0
            minParameters[1] = 1.0
        end
        if maxParameters[1] < 1.0
            minParameters[1] = 1.0
            maxParameters[1] = 1.0
        end

        # Ensure parameter ranges of constrained pairs are compatible.
        for i in 1:length(constrainedParameters)
            (pL, pR) = constrainedParameters[i]
            if minParameters[pL] > minParameters[pR]
                minParameters[pL] = minParameters[pR]
            end
            if maxParameters[pR] < maxParameters[pL]
                maxParameters[pR] = maxParameters[pL]
            end
        end

        # elite

        elite = alien(alienParameters, minParameters, maxParameters, constrainedParameters, significantFigures)
        elite.fitness = _evaluate(species, elite)

        # populationSize

        # Formula is from D. Goldberg (2002) for estimating populationSize size.
        alphabet = 2    # viz., expressions are: dominant and recessive
        schemata = significantFigures
        populationSize = Int(ceil(schemata * log(alphabet) * alphabet^schemata + log(elite.genetics.genes)))

        # generations

        generation = 1

        # Formula of D. Goldberg (2002) for estimating convergence.
        combatants = Int(max(3, populationSize÷50))
        generationsToConvergence = Int(ceil(sqrt(elite.genetics.genes)
            * log(populationSize - combatants)))

        # adults

        adults = Vector{Creature}(undef, populationSize)
        adults[1] = elite
        Threads.@threads for i in 2:populationSize
            adult = procreate(minParameters, maxParameters, constrainedParameters, significantFigures)
            adult.fitness = _evaluate(species, adult)
            adults[i] = adult
        end

        for i in 2:populationSize
            adult = adults[i]
            if adult.fitness > elite.fitness
                elite = adult
            end
        end

        # children

        children = Vector{Creature}(undef, populationSize)

        # fitness

        fitnesses = Vector{Real}(undef, populationSize)
        Threads.@threads for i in 1:populationSize
            adult = adults[i]
            fitnesses[i] = adult.fitness
        end

        new(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, constrainedParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitnesses)
    end

    function Colony(species::AbstractSpecies, probabilityOfMutation::Real, probabilityOfCrossover::Real, probabilityOfImmigrant::Real, parameterNames::Vector{String}, minParameters::Vector{Real}, maxParameters::Vector{Real}, constrainedParameters::Vector{Tuple{Int,Int}}, significantFigures::Int, populationSize::Int, generation::Int, generationsToConvergence::Int, elite::Creature, children::Vector{Creature}, adults::Vector{Creature}, fitnesses::Vector{Real})

        new(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, constrainedParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitnesses)
    end
end # Colony

# internal methods

function _tournamentPlay(c::Colony)::Creature
    combatants  = Int(max(3, c.populationSize÷50))
    contestants = Vector{Creature}(undef, combatants)
    if c.probabilityOfImmigrant > rand()
        mostFit = procreate(c.minParameters, c.maxParameters, c.constrainedParameters, c.significantFigures)
        _evaluate(c.species, mostFit)
    else
        for i in 1:combatants
            combatant = rand(1:c.populationSize)
            contestants[i] = c.adults[combatant]
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

function _mate!(c::Colony)
    # Elite creature from current generation lives into the next generation.
    c.children[1] = copy(c.elite)
    c.fitness[1]  = c.elite.fitness

    # Mate the population.
    Threads.@threads for i in 2:population
        parentA = _tournamentPlay(c)
        parentB = _tournamentPlay(c)
        child   = conceive(parentA, parentB, c.constrainedParameters, c.probabilityOfMutation, c.probabilityOfCrossover)
    end

    # Ensure their are no clones, i.e., that there are no identical twins.
    for i in 2:populationSize
        child = c.children[i]
        for j in 1:i-1
            if child == c.children[j]
                # Identical twins are not permitted within a generation.
                parentA = _tournamentPlay(c)
                parentB = _tournamentPlay(c)
                child   = conceive(parentA, parentB, c.constrainedParameters, c.probabilityOfMutation, c.probabilityOfCrossover)
                c.children[i] = child
                break
            end
        end
    end

    # Determine the fitness for each child born into the colony.
    Threads.@threads for i in 1:populationSize
        child = c.children[i]
        child.fitness  = _evaluate(c.species, child)
        c.fitnesses[i] = child.fitness
        c.children[i]  = child
    end

    return nothing
end # _mate!

function _evaluate(s::AbstractSpecies, c::Creature)::Real
    θ = parameters(c)
    modelResponse = runModel(s, θ)
    fitness = 0.0
    for i in 1:s.nExp
        sum = 0.0
        for j in 1:s.nRsp[i]
            for k in 1:s.nPts[i]
                sum = sum + ((abs(modelResponse[i][j,k] - s.resp[i][j,k])
                    / s.stdR[i][j])^θ[1])
            end
        end
        pNorm   = sum^(1.0/θ[1])
        fitness = fitness + 1.0 / pNorm
    end
    return fitness
end # _evaluate

function _stdDevElite(c::Colony)::Vector{Real}
    eliteParameters = parameters(c.elite)
    len = length(ePara)
    sde = Vector{Real}(undef, len)
    sum = zeros(Real, len)
    for i in 1:c.populationSize
        adult = c.adults[i]
        adultParameters = parameters(adult)
        for j in 1:len
            sum[j] = sum[j] + (adultParameters[j] - eliteParameters[j])^2
        end
    end
    for i in 1:len
        sde[i] = sqrt(sum[i] / (c.populationSize - 1))
    end
    return sde
end # _stdDevElite

function _2string(c::Colony, x::Real)::String
    if c.significantFigures ≤ 2
        s = @sprintf "%.1e" x;
    elseif c.significantFigures == 3
        s = @sprintf "%.2e" x;
    elseif c.significantFigures == 4
        s = @sprintf "%.3e" x;
    elseif c.significantFigures == 5
        s = @sprintf "%.4e" x;
    elseif c.significantFigures == 6
        s = @sprintf "%.5e" x;
    else
        s = @sprintf "%.6e" x;
    end
    return s
end # _2string

# exported methods

function advanceToNextGeneration!(c::Colony)
    c.generation = c.generation + 1
    _mate!(c)

    # Children become the new adults.
    Threads.@threads for i in 1:c.populationSize
        newAdult = c.children[i]
        c.adults[i] = copy(newAdult)
    end

    # Determine the elite adult for this next generation.
    c.elite = c.adults[1]
    for i in 2:populationSize
        adult = c.adult[i]
        if adult.fitness > c.elite.fitness
            c.elite = adult
        end
    end
end # advanceToNextGeneration

function report(c::Colony)::String
    s = "\n"
    s = string(s, "Statistics for generation ", c.generation)
    s = string(s, " with a population size of ", c.populationSize, ".\n")
    s = string(s, "Optimal value and population statistics for fitness:\n")
    fitness = c.elite.fitness
    if fitness ≥ 0.0
        s = string(s, "   optimal fitness  ", _2string(c, fitness), "\n")
    else
        s = string(s, "   optimal fitness ", _2string(c, fitness), "\n")
    end
    statMean = mean(c.fitness)
    if statMean ≥ 0.0
        s = string(s, "   arithmetic mean  ", _2string(c, statMean), "\n")
    else
        s = string(s, "   arithmetic mean ", _2string(c, statMean), "\n")
    end
    statMedian = median(c.fitness)
    if statMedian ≥ 0.0
        s = string(s, "   median           ", _2string(c, statMedian), "\n")
    else
        s = string(s, "   median          ", _2string(c, statMedian), "\n")
    end
    statSkewness = skewness(c.fitness)
    if statSkewness ≥ 0.0
        s = string(s, "   skewness         ", _2string(c, statSkewness), "\n")
    else
        s = string(s, "   skewness        ", _2string(c, statSkewness), "\n")
    end
    statKurtosis = kurtosis(c.fitness)
    if statKurtosis ≥ 0.0
        s = string(s, "   excess kurtosis  ", _2string(c, statKurtosis), "\n")
    else
        s = string(s, "   excess kurtosis ", _2string(c, statKurtosis), "\n")
    end
    s = string(s, "Genome from the most fit creature:\n")
    s = string(s, "   ", toString(c.elite.genetics))
    s = string(s, "Listings of [minParameter, bestParameter, maxParameter]:\n")
    bestParameters = parameters(c.elite)
    len = length(bestParameters)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", c.parameterNames[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameterNames[i], "\n")
            else
                s = string(s, i, ": ", c.parameterNames[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, "[", _2string(c.minParameters[i]), ", ",
            _2string(bestParameters[i]), ", ",
            _2string(c.maxParameters[i]), "]\n")
    end
    s = string(s, "Values for bestParameter ± error, where an RMSE ")
    s = string(s, "is computed wrt best values.\n")
    s = string(s, "Data are the parameters from all adults living ")
    s = string(s, "in the current population.\n")
    err = _stdDevElite(c)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", c.parameterNames[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameterNames[i], "\n")
            else
                s = string(s, i, ": ", c.parameterNames[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, _2string(bestParameters[i]), " ± ",
            _2string(err[i]), "\n")
    end
    s = string(s, "\n")
    return s
end # report
