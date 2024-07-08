"""
A colony is a collection of creatures in a genetic algorithm.  It is the colony
that evolves from generation to generation, leading toward an optimal condition.

A Colony is a mutable data structure as its data change from generation to
generation.  A colony has an interface of:

mutable struct Colony
    # Fields pertaining to the creatures of a colony.
    species                 # kind of creatures that comprise the colony
    probabilityOfMutation   # probability of single gene mutating at creation
    probabilityOfCrossover  # probability of gene sharing from child's parents
    probabilityOfImmigrant  # probability of an immigrant entering the colony
    parameterNames          # string representation for each model parameter
    minParameters           # lower bound when searching for optimal parameters
    maxParameters           # upper bound when searching for optimal parameters
    significantFigures      # digits of accuracy sought for optimal parameters
    # Fields pertaining to the colony itself.
    populationSize          # number of creatures comprising the colony
    generation              # counter numbering the current generation
    generationsToConvergence    # generation whose elite is the solution sought
    elite                   # most fit creature in the current generation
    children                # adults for the next generation
    adults                  # adults of the current generation
    fitness                 # fitness of each adult in the current generation
end

Constructors

    c = Colony(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, alienParameters, minParameters, maxParameters, significantFigures)

    where

        species                 An implementation of the class AbstractSpecies.
                                Here is where the user's model is introduced,
                                and where the experimental data are assigned.
        probabilityOfMutation   The probability of a gene mutation occurring at
                                conception, i.e. a gene swapping its expression.
                                Typically this is at a low value, e.g., < 0.01.
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
                                doesn't, then send θ = Vector{Float64}(undef,0).
        minParameters           Lower bounds for the parameters θ being sought.
        maxParameters           Upper bounds for the parameters θ being sought.
        significantFigures      The number of significant figures of accuracy
                                sought in a final solution for the parameters.
                                There is a strong correlation between this value
                                and how long it takes to get a solution.  Values
                                of 4 or 5 are common.  They are bound to [1, 7].

    If a minParameters[i] equals its maxParameters[i], then parameter θ[i] is
    taken to be fixed.

or the basic constructor (which typically won't be used)

    c = Colony(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitness)

IMPORTANT!!!!!

    For the parameter arrays, index 1 is reserved for selecting the best
    objective function.  Specifically, a dimensionless error is computed as
        error[i] = |model[i] - experiment[i]| / sdtDev(experiment)
    where i indexes over the set of all experimental data being considered,
    whose error is normalized by a standard deviation in the experimental
    response.  A p-norm of this experimental error is then computed, i.e.,
        ‖x‖ₚ = (∑ᵢ₌₁ⁿ|xᵢ|^p)^(1/p), p ≥ 1, with xᵢ = error[i] and n = s.nPts
    from which a quality value is tabulated as
        fitness = ∑ᵢ₌₁ⁿ 1/‖x‖ₚ[i] with n = s.nExp
    where the creature whose 'fitness' is the greatest is deemed most fit, i.e.,
    the elite.  The first entry in all parameter arrays is a 'p' for the p-norm.

    Exponent 'p', p ≥ 1, is treated here as another parameter to be optimized.
    It can be fixed or allowed to vary.  E.g., setting minParameters[1] = 1
    and maxParameters[1] = 1 would impose the taxicab norm, while setting
    minParameters[1] = 2 and maxParameters[1] = 2 would impose the Euclidean
    norm.  In addition, as 'p' goes toward infinity, the resulting p-norm will
    tend towards the infinity norm.  A setting of minParameters[1] = 1 and
    maxParameters[1] = 100 would be more typical, which would seek an optimal
    objective function (i.e., reciprocal p-norm) where the range of p lies
    within [, 100] for the specific set of data that are being considered.

Methods

advanceToNextGeneration!(c)     advances the colony to its next generation

str = report(c)                 returns a report as a string 'str' describing
                                the current generation and its statistics

"""
mutable struct Colony
    # Fields pertaining to the creatures of a colony.
    species::AbstractSpecies
    probabilityOfMutation::Float64
    probabilityOfCrossover::Float64
    probabilityOfImmigrant::Float64
    parameterNames::Vector{String}
    minParameters::Vector{Float64}
    maxParameters::Vector{Float64}
    significantFigures::Int64
    # Fields pertaining to the colony itself.
    populationSize::Int64
    generation::Int64
    generationsToConvergence::Int64
    elite::Creature
    children::Vector{Creature}
    adults::Vector{Creature}
    fitness::Vector{Float64}

    function Colony(species::AbstractSpecies, probabilityOfMutation::Float64, probabilityOfCrossover::Float64, probabilityOfImmigrant::Float64, parameterNames::Vector{String}, alienParameters::Vector{Float64}, minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64)

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
        # Their lengths should be 1 + number of parameters in the model itself.

        dim = length(parameterNames)
        str = "p in p-norm, i.e., ‖x‖ₚ = (∑ₙ₌₁ᴺ|xₙ|^p)^(1/p),  p ≥ 1"
        parameterNames[1] = str

        if (length(alienParameters) ≠ dim) && (length(alienParameters ≠ 0))
            msg = "The supplied alien has the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        if (length(minParameters) ≠ dim) || (length(maxParameters) ≠ dim)
            msg = "The supplied bounds have the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        # elite

        if length(alienParameters) > 0
            elite = alien(alienParameters, minParameters, maxParameters, significantFigures)
        else
            elite = procreate(minParameters, maxParameters, significantFigures)
        end
        elite.fitness = _evaluate(species, elite)

        # populationSize

        # Formula is from D. Goldberg (2002) for estimating populationSize size.
        alphabet = 2    # viz., expressions dominant and recessive
        schemata = significantFigures
        populationSize = Int64(ceil(schemata * log(alphabet) * alphabet^schemata + log(elite.genetics.genes)))

        # generations

        generation = 1

        # Formula of D. Goldberg (2002) for estimating convergence.
        combatants = Int64(max(3, populationSize÷50))
        generationsToConvergence = Int64(ceil(sqrt(elite.genetics.genes)
            * log(populationSize - combatants)))

        # adults

        adults = Vector{Creature}(undef, populationSize)
        adults[1] = elite
        Threads.@threads for i in 2:populationSize
            adult = procreate(minParameters, maxParameters, significantFigures)
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

        fitness = Vector{Float64}(undef, populationSize)
        Threads.@threads for i in 1:populationSize
            adult = adults[i]
            fitness[i] = adult.fitness
        end

        new(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitness)
    end

    function Colony(species::AbstractSpecies, probabilityOfMutation::Float64, probabilityOfCrossover::Float64, probabilityOfImmigrant::Float64, parameterNames::Vector{String}, minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64, populationSize::Int64, generation::Int64, generationsToConvergence::Int64, elite::Creature, children::Vector{Creature}, adults::Vector{Creature}, fitness::Vector{Float64})

        new(species, probabilityOfMutation, probabilityOfCrossover, probabilityOfImmigrant, parameterNames, minParameters, maxParameters, significantFigures, populationSize, generation, generationsToConvergence, elite, children, adults, fitness)
    end
end # Colony

# internal methods

function _tournamentPlay(c::Colony)::Creature
    combatants  = Int64(max(3, c.populationSize÷50))
    contestants = Vector{Creature}(undef, combatants)
    if c.probabilityOfImmigrant > rand()
        mostFit = Creature(c.minParameters, c.maxParameters, c.significantFigures)
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
        child   = conceive(parentA, parentB, probabilityOfMutation, probabilityOfCrossover)
    end

    # Ensure their are no clones, i.e., that there are no identical twins.
    for i in 2:populationSize
        child = c.children[i]
        for j in 1:i-1
            if child == c.children[j]
                # Identical twins are not permitted within a generation.
                parentA = _tournamentPlay(c)
                parentB = _tournamentPlay(c)
                c.children[i] = conceive(parentA, parentB, probabilityOfMutation, probabilityOfCrossover)
                j = i - 1
            end
        end
    end

    # Determine the fitness for each child born into the colony.
    Threads.@threads for i in 1:populationSize
        child = c.children[i]
        child.fitness = _evaluate(c.species, child)
        c.fitness[i]  = child.fitness
        c.children[i] = child
    end

    return nothing
end # _mate!

function _evaluate(s::AbstractSpecies, c::Creature)::Float64
    θ = parameters(c)
    # Ensure the p-norm is bound from below such that p ≥ 1.
    if θ[1] ≥ 1.0
        p = θ[1]
    else
        p = 1.0
    end
    modR = runModel(s, θ)
    fitness = 0.0
    for i in 1:s.nExp
        sum = 0.0
        for j in 1:s.nRsp[i]
            for k in 1:s.nPts[i]
                sum = sum + ((abs(modR[i][j,k] - s.resp[i][j,k])
                    / s.stdR[i][j])^p)
            end
        end
        pNorm   = sum^(1.0/p)
        fitness = fitness + 1.0 / pNorm
    end
    return fitness
end # _evaluate

function _stdDevElite(c::Colony)::Vector{Float64}
    bst = parameters(c.elite)
    len = length(bst)
    sde = Vector{Float64}(undef, len)
    sum = zeros(Float64, len)
    for i in 1:c.populationSize
        adult = c.adults[i]
        par_i = parameters(adult)
        for j in 1:len
            sum[j] = sum[j] + (par_i[j] - bst[j])^2
        end
    end
    for i in 1:len
        sde[i] = sqrt(sum[i] / (c.populationSize - 1))
    end
    return sde
end # _stdDevElite

function _2string(c::Colony, x::Float64)::String
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
    statSkew = skewness(c.fitness)
    if statSkew ≥ 0.0
        s = string(s, "   skewness         ", _2string(c, statSkew), "\n")
    else
        s = string(s, "   skewness        ", _2string(c, statSkew), "\n")
    end
    statKurtosis = kurtosis(c.fitness)
    if statKurtosis ≥ 0.0
        s = string(s, "   excess kurtosis  ", _2string(c, statKurtosis), "\n")
    else
        s = string(s, "   excess kurtosis ", _2string(c, statKurtosis), "\n")
    end
    s = string(s, "Genetic genome from the most fit creature:\n")
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
