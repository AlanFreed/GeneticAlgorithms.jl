"""
A colony is a collection of creatures in a genetic algorithm.  It is the colony
that evolves from generation to generation, leading toward an optimum condition.

A colony has an interface of:

struct Colony{Parameters} where Parameters <: AbstractParameters
    parameters              # model parameters: object from AbstractParameters
    data                    # experimental data that model is to be fit against
    # Fields pertaining to the creatures of a colony.
    probability_mutation    # probability of single gene mutating at creation
    probability_crossover   # probability of a child sharing its parents genes
    probability_immigrant   # probability of an immigrant entering the colony
    parameters_names        # string representation for each model parameter
    parameters_min          # lower bounds when searching for optimal parameters
    parameters_max          # upper bounds when searching for optimal parameters
    parameters_constrained  # tuples specifying constraints of: left < right
    significant_figures     # digits of accuracy sought for optimal parameters
    # Fields pertaining to the colony itself.
    population_size         # number of creatures comprising the colony
    generation              # counter specifying the current generation
    generations_to_convergence  # generation whose elite is the solution sought
    elite                   # most fit resident in the current generation
    children                # adults for the next generation
    adults                  # adults from the current generation
    fitness                 # fitness of each adult in the current generation
end

Constructors

    The constructor most likely to be called by a user.

    c = Colony(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_alien, parameters_min, parameters_max, parameters_constrained, significant_figures)

    where

        parameters              An object holding the model's parameters as a
                                𝑚𝑢𝑡𝑎𝑏𝑙𝑒 composite type whose parent type is
                                AbstractParameters.  What is important here is
                                the object, not the values held by the object.
        data                    Experimental data that the model is to be fit
                                against, whose optimal parameters θ are sought.
        probability_mutation    The probability of a gene mutation occurring at
                                conception, i.e. a gene swapping its expression.
                                Typically this is a low value, e.g., < 0.01.
        probability_crossover   The probability of a crossover event occurring
                                at conception, i.e., a chromosome splitting.
                                Typically this is a large value, e.g., > 0.85.
        probability_immigrant   The probability of introducing an immigrant
                                into the gene pool after the first generation.
                                (The first generation is all immigrants.)  Use
                                as you see fit.  I usually assign it so as to
                                introduce one immigrant per generation or so.
        parameters_alien        If the user has a best guess for θ, then this is
                                where these parameters are input.  If the user
                                doesn't, then send: θ = Vector{Real}(undef,0).
        parameters_min          Lower bounds for the parameters θ being sought.
        parameters_max          Upper bounds for the parameters θ being sought.
        parameters_constrained  Tuples of indices (left, right) that impose an
                                inequality constraint, viz.: θ[left] < θ[right].
                                No constraints? Send: Vector{Tuples}(undef,0).
        significant_figures     The number of significant figures of accuracy
                                sought in a final solution for the parameters.
                                There is a strong correlation between this value
                                and how long it takes to get a solution.  Values
                                of 4 or 5 are common.  They are bound to [1, 7].

    If parameter parameters_min[i] equals parameter parameters_max[i] for any
    index i, then parameter θ[i] is taken to be fixed.

Methods

advance_to_next_generation!(c)  Advances the colony 'c' to its next generation.

str = report(c)                 Returns a human-readable report via a string
                                'str,' which describes health of the current
                                generation and its statistics.
"""
struct Colony
    parameters::AbstractParameters
    data::ExperimentalData
    # Fields pertaining to the creatures of a colony.
    probability_mutation::Real
    probability_crossover::Real
    probability_immigrant::Real
    parameters_name::Vector{String}
    parameters_min::Vector{Real}
    parameters_max::Vector{Real}
    parameters_constrained::Vector{Tuple{Integer,Integer}}
    significant_figures::Integer
    # Fields pertaining to the colony itself.
    population_size::Integer
    generations_to_convergence::Integer
    # Mutable fields pertaining to the colony itself.
    generation::Counter
    elite::Creature
    children::Vector{Creature}
    adults::Vector{Creature}
    fitness::Vector{Real}

    # constructor

    function Colony(parameters::AbstractParameters, data::ExperimentalData, probability_mutation::Real, probability_crossover::Real, probability_immigrant::Real, parameters_alien::Vector{Real}, parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer)

        # bound inputs

        if significant_figures < 1
            significant_figures = 1
        end
        if significant_figures > 7
            significant_figures = 7
        end

        if probability_mutation < 1.0e-6
            probability_mutation = 1.0e-6
        end
        if probability_mutation > 0.999999
            probability_mutation = 0.999999
        end

        if probability_crossover < 1.0e-6
            probability_crossover = 1.0e-6
        end
        if probability_crossover > 0.999999
            probability_crossover = 0.999999
        end

        if probability_immigrant < 1.0e-6
            probability_immigrant = 1.0e-6
        end
        if probability_immigrant > 0.999999
            probability_immigrant = 0.999999
        end

        # Verify dimensions.

        N = fieldcount(typeof(parameters))

        if (length(parameters_alien) ≠ N) && (length(parameters_alien) ≠ 0)
            msg = "The supplied alien has the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        if (length(parameters_min) ≠ N) || (length(parameters_max) ≠ N)
            msg = "The supplied bounds have the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        # Ensure the parameters range of constrained pairs are compatible.

        for i in 1:length(parameters_constrained)
            (pL, pR) = parameters_constrained[i]
            if parameters_min[pL] > parameters_min[pR]
                parameters_min[pL] = parameters_min[pR]
            end
            if parameters_max[pR] < parameters_max[pL]
                parameters_max[pR] = parameters_max[pL]
            end
        end

        #=
        # Verify that the fields of object parameters are all Real valued.

        for n in 1:N
            symbol  = fieldname(typeof(parameters), n)
            if !isa(fieldtype(typeof(parameters), symbol), Real)
                msg = "All fields in object parameters must belong to Real."
                error(msg)
            end
        end
        =#

        # Get parameter names from the fields of object parameters.

        parameters_name = Vector{String}(undef, N)
        for n in 1:N
            symbol = fieldname(typeof(parameters), n)
            parameters_name[n] = String(symbol)
        end

        # Create the model.

        model = Model(parameters, data)

        # elite

        elite = alien(parameters_alien, parameters_min, parameters_max, parameters_constrained, significant_figures)
        set!(model, phenotypes(elite))
        set!(elite.fitness, _evaluate(model))

        # population_size

        # Formula is from D. Goldberg (2002) for estimating population size.
        alphabet = 2    # viz., expressions are: dominant and recessive
        schemata = significant_figures
        population_size = Int(ceil(schemata * log(alphabet) * alphabet^schemata
            + log(elite.genetics.genes)))

        # generations

        generation = Counter(1)

        # Formula is from D. Goldberg (2002) for estimating convergence.
        combatants = Int(max(3, population_size÷50))
        generations_to_convergence = Int(ceil(sqrt(elite.genetics.genes)
            * log(population_size - combatants)))

        # adults

        adults = Vector{Creature}(undef, population_size)
        adults[1] = elite
        Threads.@threads for i in 2:population_size
            adult = procreate(parameters_min, parameters_max, parameters_constrained, significant_figures)
            thread_model = Model(parameters, data)
            set!(thread_model, phenotypes(adult))
            set!(adult.fitness, _evaluate(thread_model))
            adults[i] = adult
        end

        for i in 2:population_size
            adult = adults[i]
            if adult.fitness > elite.fitness
                elite = adult
            end
        end

        # children

        children = Vector{Creature}(undef, population_size)

        # fitness

        fitness = Vector{Real}(undef, population_size)
        for i in 1:population_size
            adult = adults[i]
            fitness[i] = get(adult.fitness)
        end

        new(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_name, parameters_min, parameters_max, parameters_constrained, significant_figures, population_size, generations_to_convergence, generation, elite, children, adults, fitness)
    end

    function Colony(parameters::AbstractParameters, data::ExperimentalData, probability_mutation::Real, probability_crossover::Real, probability_immigrant::Real, parameters_name::Vector{String}, parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer, population_size::Integer, generations_to_convergence::Integer, generation::Counter, elite::Creature, children::Vector{Creature}, adults::Vector{Creature}, fitness::Vector{Real})

        new(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_name, parameters_min, parameters_max, parameters_constrained, significant_figures, population_size, generations_to_convergence, generation, elite, children, adults, fitness)
    end
end # Colony

# internal methods

function _tournamentPlay(c::Colony)::Creature
    combatants  = Int(max(3, c.population_size÷50))
    contestants = Vector{Creature}(undef, combatants)
    if c.probability_immigrant > rand()
        mostfit = procreate(c.parameters_min, c.parameters_max, c.parameters_constrained, c.significant_figures)
        model = Model(c.parameters, c.data)
        set!(model, phenotypes(mostfit))
        set!(mostfit.fitness, _evaluate(model))
    else
        for i in 1:combatants
            combatant = rand(1:c.population_size)
            contestants[i] = c.adults[combatant]
        end
        mostfit = contestants[1]
        for i in 2:combatants
            contestant = contestants[i]
            if contestant.fitness > mostfit.fitness
                mostfit = contestant
            end
        end
    end
    return mostfit
end # _tournamentPlay

function _mate!(c::Colony)
    # Elite creature from current generation lives into the next generation.
    c.children[1] = copy(c.elite)
    c.fitness[1]  = get(c.elite.fitness)

    # Mate the population.
    Threads.@threads for i in 2:c.population_size
        parentA = _tournamentPlay(c)
        parentB = _tournamentPlay(c)
        count = 1
        while parentB == parentA
            count   = count + 1
            parentB = _tournamentPlay(c)
            if count == 25
                # This is for safety. It should not occur in practice.
                break
            end
        end

        child   = conceive(parentA, parentB, c.parameters_constrained, c.probability_mutation, c.probability_crossover)
    end

    # Ensure their are no clones, i.e., that there are no identical twins.
    for i in 2:population_size
        child = c.children[i]
        for j in 1:i-1
            if child == c.children[j]
                # Identical twins are not permitted within a generation.
                parentA = _tournamentPlay(c)
                parentB = _tournamentPlay(c)
                child   = conceive(parentA, parentB, c.parameters_constrained, c.probability_mutation, c.probability_crossover)
                c.children[i] = child
                break
            end
        end
    end

    # Determine the fitness for each child born into the colony.
    Threads.@threads for i in 2:population_size
        child = c.children[i]
        model = Model(c.parameters, c.data)
        set!(model, phenotypes(child))
        set!(child.fitness, _evaluate(model))
        c.fitness[i]  = get(child.fitness)
        c.children[i] = child
    end

    return nothing
end # _mate!

function _evaluate(m::Model)::Real
    model_responses = solve(m)
    fitness = 0.0
    for i in 1:m.d.experiments
        squared_error = 0.0
        for j in 1:m.d.variables_response[i]
            for k in 1:m.d.data_points[i]
                squared_error = (squared_error
                    + ((model_responses[i][j,k] - m.d.responses[i][j,k])
                    / m.d.responses_std[i][j])^2)
            end
        end
        twonorm = sqrt(squared_error)
        fitness = fitness + 1.0 / twonorm
    end
    return fitness
end # _evaluate

function _stdDevElite(c::Colony)::Vector{Real}
    parameters_elite = phenotypes(c.elite)
    len = length(parameters_elite)
    sde = Vector{Real}(undef, len)
    sum = zeros(Real, len)
    for i in 1:c.population_size
        adult = c.adults[i]
        parameters_adult = phenotypes(adult)
        for j in 1:len
            sum[j] = sum[j] + (parameters_adult[j] - parameters_elite[j])^2
        end
    end
    for i in 1:len
        sde[i] = sqrt(sum[i] / (c.population_size - 1))
    end
    return sde
end # _stdDevElite

function _2string(c::Colony, x::Real)::String
    if c.significant_figures ≤ 2
        s = @sprintf "%.1e" x;
    elseif c.significant_figures == 3
        s = @sprintf "%.2e" x;
    elseif c.significant_figures == 4
        s = @sprintf "%.3e" x;
    elseif c.significant_figures == 5
        s = @sprintf "%.4e" x;
    elseif c.significant_figures == 6
        s = @sprintf "%.5e" x;
    else
        s = @sprintf "%.6e" x;
    end
    return s
end # _2string

# exported methods

function advance_to_next_generation!(c::Colony)
    set!(c.generation, get(c.generation)+1)
    _mate!(c)

    # Children become the new adults.
    for i in 1:c.population_size
        c.adults[i] = copy(c.children[i])
    end

    # Determine the elite adult for this next generation.
    elite = c.adults[1]
    for i in 2:population_size
        adult = c.adult[i]
        if adult.fitness > c.elite.fitness
            elite = adult
        end
    end
    set!(c.elite.fitness, get(elite.fitness))
    for chromosome in 1:c.elite.genetics.chromosomes
        genotype_c = c.elite.genetics.genotypes[chromosome]
        genotype_e =   elite.genetics.genotypes[chromosome]
        for gene in 1:genotype_c.genes
            setindex!(genotype_c, getindex(genotype_e, gene), gene)
        end
    end
    return nothing
end # advanceToNextGeneration

function report(c::Colony)::String
    s = "\n"
    s = string(s, "Statistics for generation ", c.generation)
    s = string(s, " with a population size of ", c.population_size, ".\n")
    s = string(s, "Optimum fitness and population statistics for fitness:\n")
    fitness = get(c.elite.fitness)
    if fitness ≥ 0.0
        s = string(s, "   optimum fitness  ", _2string(c, fitness), "\n")
    else
        s = string(s, "   optimum fitness ", _2string(c, fitness), "\n")
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
    s = string(s, "The genome from the most fit creature:\n")
    s = string(s, "   ", tostring(c.elite.genetics))
    s = string(s, "Lists for [parameter_min, parameter_best, parameter_max]:\n")
    parameters_best = phenotypes(c.elite)
    len = length(parameters_best)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", c.parameters_name[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameters_name[i], "\n")
            else
                s = string(s, i, ": ", c.parameters_name[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, "[", _2string(c, c.parameters_min[i]), ", ",
            _2string(c, parameters_best[i]), ", ",
            _2string(c, c.parameters_max[i]), "]\n")
    end
    s = string(s, "Values for best parameter ± error, where an RMSE ")
    s = string(s, "is computed wrt best values.\n")
    s = string(s, "Data are the parameters from all adults living ")
    s = string(s, "in the current population.\n")
    err = _stdDevElite(c)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", c.parameters_name[i], "\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameters_name[i], "\n")
            else
                s = string(s, i, ": ", c.parameters_name[i], "\n")
            end
            s = string(s, "       ")
        end
        s = string(s, _2string(c, parameters_best[i]), " ± ",
            _2string(c, err[i]), "\n")
    end
    s = string(s, "\n")
    return s
end # report
