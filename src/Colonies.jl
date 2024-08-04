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
    use_fitness             # select which fitness type to use, ∈ [1, 4]
end

Constructors

    The constructor most likely to be called by a user.

    c = Colony(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_alien, parameters_min, parameters_max, parameters_constrained, significant_figures, use_fitness)

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
        use_fitness             This selects the objective function to use from
                                four possible choices with the default set at 3:
                                    1) minimizes expectation in absolute error
                                    2) minimizes expectation in squared error
                                    3) minimizes variance in error
                                    4) maximizes covariance between experiment
                                        and model

    If parameter parameters_min[i] equals parameter parameters_max[i] for any
    index i, then parameter θ[i] is taken to be fixed.

Methods

advance_to_next_generation!(c)  Advances colony 'c' to its next generation.

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
    fitness::Vector{Vector{Real}}
    use_fitness::Integer

    # constructor

    function Colony(parameters::AbstractParameters, data::ExperimentalData, probability_mutation::Real, probability_crossover::Real, probability_immigrant::Real, parameters_alien::Vector{Real}, parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer, use_fitness::Integer = 3)

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

        # Verify that the fields of object parameters are all Real valued.

        for n in 1:N
            symbol  = fieldname(typeof(parameters), n)
            if !(fieldtype(typeof(parameters), symbol) isa Real)
                msg = "All fields in object parameters must belong to Real."
                error(msg)
            end
        end
        =#

        # Get parameter names from the fields of object parameters.

        parameters_name = Vector{String}(undef, N)
        for n in 1:N
            symbol = fieldname(typeof(parameters), n)
            parameters_name[n] = chomp(String(symbol))
        end

        # Create the model.

        model = Model(parameters, data)

        # elite

        elite = alien(parameters_alien, parameters_min, parameters_max, parameters_constrained, significant_figures)
        set!(model, phenotypes(elite))
        fitness = _evaluate(model)
        for i in 1:fitness_types
            elite.fitness[i] = fitness[i]
        end

        # population_size

        # Formula is from D. Goldberg (2002) for estimating population size.
        alphabet = 2    # viz., expressions are: dominant and recessive
        schemata = significant_figures
        population_size = Int(ceil(alphabet^schemata * schemata * log(alphabet) 
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
            thread_fitness = _evaluate(thread_model)
            for j in 1:fitness_types
                adult.fitness[j] = thread_fitness[j]
            end
            adults[i] = adult
        end

        for i in 2:population_size
            adult = adults[i]
            if adult.fitness[use_fitness] > elite.fitness[use_fitness]
                elite = adult
            end
        end

        # children

        children = Vector{Creature}(undef, population_size)

        # fitness

        fitness = Vector{Vector}(undef, population_size)

        for i in 1:population_size
            adult = adults[i]
            fitness_adult = Vector{Real}(undef, fitness_types)
            for j in 1:fitness_types
                fitness_adult[j] = adult.fitness[j]
            end
            fitness[i] = fitness_adult
        end

        new(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_name, parameters_min, parameters_max, parameters_constrained, significant_figures, population_size, generations_to_convergence, generation, elite, children, adults, fitness, use_fitness)
    end

    function Colony(parameters::AbstractParameters, data::ExperimentalData, probability_mutation::Real, probability_crossover::Real, probability_immigrant::Real, parameters_name::Vector{String}, parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_constrained::Vector{Tuple{Integer,Integer}}, significant_figures::Integer, population_size::Integer, generations_to_convergence::Integer, generation::Counter, elite::Creature, children::Vector{Creature}, adults::Vector{Creature}, fitness::Vector{Vector{Real}}, use_fitness::Integer = 3)

        new(parameters, data, probability_mutation, probability_crossover, probability_immigrant, parameters_name, parameters_min, parameters_max, parameters_constrained, significant_figures, population_size, generations_to_convergence, generation, elite, children, adults, fitness, use_fitness)
    end
end # Colony

# internal methods

function _tournamentPlay(c::Colony)::Creature
    combatants  = Int(max(3, c.population_size÷50))
    contestants = Vector{Creature}(undef, combatants)
    if c.probability_immigrant > rand()
        # Add an immigrant to the population.
        mostfit = procreate(c.parameters_min, c.parameters_max, c.parameters_constrained, c.significant_figures)
        model = Model(c.parameters, c.data)
        set!(model, phenotypes(mostfit))
        fitness = _evaluate(model)
        for i in 1:fitness_types
            mostfit.fitness[i] = fitness[i]
        end
    else
        # Select a creature for mating.
        for i in 1:combatants
            combatant = rand(1:c.population_size)
            contestants[i] = c.adults[combatant]
        end
        mostfit = contestants[1]
        for i in 2:combatants
            contestant = contestants[i]
            if (contestant.fitness[c.use_fitness]
                > mostfit.fitness[c.use_fitness])
                mostfit = contestant
            end
        end
    end
    return mostfit
end # _tournamentPlay

function _mate!(c::Colony)
    # Elite creature from current generation lives into the next generation.
    c.children[1] = copy(c.elite)

    # Mate the population.
    Threads.@threads for i in 2:c.population_size
        parentA = _tournamentPlay(c)
        parentB = _tournamentPlay(c)
        self_conception = 0
        while parentB == parentA
            self_conception = self_conception + 1
            parentB = _tournamentPlay(c)
            if self_conception == 25
                # This is for safety. It should not occur in practice.
                break
            end
        end
        child = conceive(parentA, parentB, c.parameters_constrained, c.probability_mutation, c.probability_crossover)
        c.children[i] = child
    end

    # Ensure there are no clones, i.e., that there are no identical twins.
    for i in 2:c.population_size
        child = c.children[i]
        for j in 1:i-1
            if child == c.children[j]
                # Identical twins are not permitted within a generation.
                parentA = _tournamentPlay(c)
                parentB = _tournamentPlay(c)
                new_child = conceive(parentA, parentB, c.parameters_constrained, c.probability_mutation, c.probability_crossover)
                c.children[i] = new_child
                break
            end
        end
    end

    # Determine the fitness for each child born into the colony.
    Threads.@threads for i in 2:c.population_size
        child = c.children[i]
        model = Model(c.parameters, c.data)
        set!(model, phenotypes(child))
        fitness = _evaluate(model)
        for j in 1:fitness_types
            child.fitness[j] = fitness[j]
            c.fitness[i][j]  = fitness[j]
        end
        c.children[i] = child
    end

    return nothing
end # _mate!

# Statistical functions used to evaluate the fitness of a creature.

# Statistical moments

function _first_moment(v::Vector{Real})::Real
    N = length(v)
    sum = 0.0
    for n in 1:N
        sum = sum + v[n]
    end
    moment = sum / N
    return moment
end # _first_moment

function _second_moment(v::Vector{Real})::Real
    N = length(v)
    sum = 0.0
    for n in 1:N
        sum = sum + v[n]^2
    end
    moment = sum / N
    return moment
end # _second_moment

function _mixed_moment(v1::Vector{Real}, v2::Vector{Real})::Real
    N = length(v1)
    sum = 0.0
    for n in 1:N
        sum = sum + v1[n] * v2[n]
    end
    moment = sum / N
    return moment
end # _mixed_moment

# Sample variance

function _sample_variance(v::Vector{Real})::Real
    variance = _second_moment(v) - _first_moment(v)^2
    return variance
end # _sample_variance

# Sample covariance

function _sample_covariance(v1::Vector{Real}, v2::Vector{Real})::Real
    covariance = _mixed_moment(v1, v2) - _first_moment(v1) * _first_moment(v2)
    return covariance
end # _sample_covariance

# Objective functions

function _objfn(response_experiment::Vector{Real}, response_model::Vector{Real})::Vector{Real}

    if length(response_experiment) == length(response_model)
        N = length(response_experiment)
    else
        msg = "Lengths for model and experimental response vectors must equal."
        error(msg)
    end

    # Determine the factor of normalization.
    res_max = 0.0
    for n in 1:N
        res_max = max(res_max, abs(response_experiment[n]))
    end

    # Normalize the incoming vectors.
    res_exp = Vector{Real}(undef, N)
    res_mod = Vector{Real}(undef, N)
    res_err = Vector{Real}(undef, N)
    for n in 1:N
        res_exp[n] = response_experiment[n] / res_max   # ∈ [-1, 1]
        res_mod[n] = response_model[n] / res_max
        res_err[n] = abs(res_exp[n] - res_mod[n])
    end

    # Get the various statistics for these normalized vectors.
    E_exp   = _first_moment(res_exp)
    E_mod   = _first_moment(res_mod)
    VAR_exp = _sample_variance(res_exp)
    VAR_mod = _sample_variance(res_mod)
    COV     = _sample_covariance(res_exp, res_mod)

    # (minimize) expectation for the magnitude of error
    ϕ₁ = _first_moment(res_err)

    # Here error is the experimental response minus the model's response.

    # (minimize) expectation for the squared error:
    # E((res_exp - res_mod)²)
    ϕ₂ = VAR_exp + VAR_mod - 2COV + E_exp^2 + E_mod^2 -2E_exp * E_mod

    # (minimize) variance in the error:
    # VAR(res_exp - res_mod) where VAR(X) = E(X²) - (E(X))²
    ϕ₃ = VAR_exp + VAR_mod - 2COV

    # (maximize) covariance between experiment response and model prediction
    ϕ₄ = COV

    # Reciprocal values will maximize minimum objective functions.
    fitness = Vector{Real}(undef, fitness_types)
    fitness[1] = 1 / ϕ₁
    fitness[2] = 1 / ϕ₂
    fitness[3] = 1 / ϕ₃
    fitness[4] = ϕ₄

    return fitness
end # _objfn

function _evaluate(m::Model)::Vector{Real}

    model_responses = solve(m)

    data_points = 0
    for exp in 1:m.d.experiments
        data_points = (data_points
            + m.d.variables_response[exp] * m.d.data_points[exp])
    end

    fitness = Vector{Real}(undef, fitness_types)
    for i in 1:m.d.experiments
        for j in 1:m.d.variables_response[i]
            ϕ = _objfn(m.d.responses[i][j,:], model_responses[i][j,:])
            ratio = m.d.data_points[i] / data_points
            for k in 1:fitness_types
                fitness[k] = fitness[k] + ratio * ϕ[k]
            end
        end
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
    return chomp(s)
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
    for i in 2:c.population_size
        adult = c.adults[i]
        if adult.fitness[c.use_fitness] > elite.fitness[c.use_fitness]
            elite = adult
        end
    end
    for i in 1:fitness_types
        c.elite.fitness[i] = elite.fitness[i]
    end
    for chromosome in 1:c.elite.genetics.chromosomes
        (c.elite.genetics.genotypes[chromosome]
            = elite.genetics.genotypes[chromosome])
    end
    return nothing
end # advanceToNextGeneration

function report(c::Colony)::String
    s = ""
    s = string(s, "Statistics for generation ", tostring(c.generation))
    s = string(s, " with a population size of ", c.population_size, ".\n")
    s = string(s, "Optimum fitness and population statistics for fitness,\n")
    if c.use_fitness < 4
        s = string(s, "where fitness minimizes ")
    else
        s = string(s, "where fitness maximizes ")
    end
    if c.use_fitness == 1
        s = string(s, "expectation for the magnitude in error, are:\n")
    elseif c.use_fitness == 2
        s = string(s, "expectation for the squared error, are:\n")
    elseif c.use_fitness == 3
        s = string(s, "sample variance for the error, are:\n")
    else # c.use_fitness == 4
        s = string(s, "sample covariance between experiment and model, are:\n")
    fitness = c.elite.fitness[c.use_fitness]
    if fitness ≥ 0.0
        s = string(s, "   optimum fitness  ", _2string(c, fitness), "\n")
    else
        s = string(s, "   optimum fitness ", _2string(c, fitness), "\n")
    end
    stat_mean = mean(c.fitness[:][c.use_fitness])
    if stat_mean ≥ 0.0
        s = string(s, "   arithmetic mean  ", _2string(c, stat_mean), "\n")
    else
        s = string(s, "   arithmetic mean ", _2string(c, stat_mean), "\n")
    end
    stat_median = median(c.fitness[:][c.use_fitness])
    if stat_median ≥ 0.0
        s = string(s, "   median           ", _2string(c, stat_median), "\n")
    else
        s = string(s, "   median          ", _2string(c, stat_median), "\n")
    end
    stat_skewness = skewness(c.fitness[:][c.use_fitness])
    if stat_skewness ≥ 0.0
        s = string(s, "   skewness         ", _2string(c, stat_skewness), "\n")
    else
        s = string(s, "   skewness        ", _2string(c, stat_skewness), "\n")
    end
    stat_kurtosis = kurtosis(c.fitness[:]{c.use_fitness})
    if stat_kurtosis ≥ 0.0
        s = string(s, "   excess kurtosis  ", _2string(c, stat_kurtosis), "\n")
    else
        s = string(s, "   excess kurtosis ", _2string(c, stat_kurtosis), "\n")
    end
    s = string(s, "The genome from the most fit creature:\n")
    s = string(s, tostring(c.elite.genetics), "\n")
    s = string(s, "Lists for [parameter_min, parameter_best, parameter_max]:\n")
    parameters_best = phenotypes(c.elite)
    len = length(parameters_best)
    for i in 1:len
        if len < 10
            s = string(s, i, ": ", c.parameters_name[i], " ∈\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameters_name[i], " ∈\n")
            else
                s = string(s, i, ": ", c.parameters_name[i], " ∈\n")
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
            s = string(s, i, ": ", c.parameters_name[i], " =\n")
            s = string(s, "      ")
        else
            if i < 10
                s = string(s, i, ":  ", c.parameters_name[i], " =\n")
            else
                s = string(s, i, ": ", c.parameters_name[i], " =\n")
            end
            s = string(s, "       ")
        end
        s = string(s, _2string(c, parameters_best[i]), " ± ",
            _2string(c, err[i]), "\n")
    end
    s = string(s, "\n")
    return s
end # report
