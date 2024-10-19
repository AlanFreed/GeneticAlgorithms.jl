"""
A colony is a population of creatures in a genetic algorithm. It is the colony that evolves from one generation to the next, converging toward an optimal state of fitness.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# Colony

```julia
struct Colony
    myparameters::AbstractParameters
    mydata::ExperimentalData
    # Fields pertaining to the creatures of a colony.
    probability_mutation::Float64
    probability_crossover::Float64
    probability_immigrant::Float64
    parameters_min::Vector{PF.PhysicalScalar}
    parameters_max::Vector{PF.PhysicalScalar}
    parameters_constrained::Vector{Tuple{Int,Int}}
    significant_figures::Int
    # Fields pertaining to the colony itself.
    population::Int
    generations::Int
    elite::Creature
    # Mutable fields pertaining to the colony itself.
    generation::PF.MInteger
    children::Vector{Creature}
    adults::Vector{Creature}
end
```

## Constructors

The constructor most likely to be called by a user.

```julia
colony = Colony(myparameters::AbstractParameters,
                mydata::TheData, 
                # Fields pertaining to the creatures of a colony.
                probability_mutation::Float64, 
                probability_crossover::Float64, 
                probability_immigrant::Float64,
                parameters_alien::Vector{PF.PhysicalScalar},
                parameters_min::Vector{PF.PhysicalScalar},
                parameters_max::Vector{PF.PhysicalScalar}, 
                parameters_constrained::Vector{Tuple{Int,Int}}, 
                significant_figures::Int)
```
where
1. *myparameters* is an object holding the model's parameters θ whose supertype is AbstractParameters. What is needed by this constructor is the object itself, i.e., its fields for meta programming, not the actual values held by these fields.
2. *mydata* holds the experimental data that the model is to be fit against, whose optimal parameters θ are sought, and is a container for the model's predictions, too.
3. *probability_mutation* is the probability of a gene mutation occurring at conception, i.e., a gene swapping its expression. Typically this is at a low value, e.g., < 0.01.
4. *probability_crossover* is the probability of a crossover event occurring at conception, i.e., the splitting and pairing of chromosomes coming from two parents. Typically this is at a large value, e.g., > 0.99.
5. *probability_immigrant* is the probability of introducing an immigrant into tournament play (mate selection) after the first generation.
6. *parameters_alien* is the user's best guess for parameters θ. If the user doesn't have a best guess, then the user should send `parameters_alien = Vector{PF.PhysicalScalar}(undef, 0)` to the constructor.
7. *parameters_min* are lower bounds for the parameters θ being sought.
8. *parameters_max* are upper bounds for the parameters θ being sought.
9. *parameters_constrained* are tuples of indices (left, right) that impose an inequality constraint, viz.: `θ[left] < θ[right]`. If there are no constraints, then the user should send `parameters_constrained = Vector{Tuples}(undef, 0)` to the constructor.
10. *significant_figures* is the number of significant figures of accuracy sought in a final solution for the parameters. There is a strong correlation between this value and how long it takes to get a solution. Values of 4 or 5 are common. They are bound to the interval [1…8].

> If a minimum parameter equals its maximum parameter at any given index in the sense that `parameters_min[i] ≈ parameters_max[i]` at say index i, then that model parameter, i.e., θ[i], will be taken to be fixed at its specified value.

### Full constructor (used by JSON3)

```julia
colony = Colony(myparameters::AbstractParameters,
                mydata::TheData,
                # Fields pertaining to the creatures of a colony.
                probability_mutation::Float64,
                probability_crossover::Float64,
                probability_immigrant::Float64,
                parameters_min::Vector{PF.PhysicalScalar},
                parameters_max::Vector{PF.PhysicalScalar},
                parameters_constrained::Vector{Tuple{Int,Int}},
                significant_figures::Int,
                # Fields pertaining to the colony itself.
                population::Int,
                generations::Int,
                elite::Creature,
                # Mutable fields pertaining to the colony itself.
                generation::PF.MInteger,
                children::Vector{Creature},
                adults::Vector{Creature})
```
whose arguments, in addition to those of the first constructor, are:
1. *population* establishes how many creatures comprise the colony. This is fixed from one generation to the next, and is established through a formula derived from probability theory by D. Goldberg (2002) in the other constructor.
2. *generations* specifies the number of generations estimated to achieve a convergence for the parameters sought, in a statistical sense. This is established through a formula derived from probability theory by D. Goldberg (2002) in the other constructor.
3. *elite* is the most fit creature in the population at its current generation.
4. *generation* is a counter representing the current generation along a solution path that will end when `generation == generations`.
5. *children* is a vector holding the newly conceived creatures of a generation that are to become the adults of the next generation.
6. *adults* is a vector of creatures comprising the current generation, and who can mate.

## Methods

### Advance the colony of creatures to its next generation.

```julia
advance!(colony::Colony)
```
> increments a solution along its path, advancing it by one generation.

### Write a report

```julia
text = report(colony::Colony)
```
> returns a human-readable report via a string `text` that describes the health of the current generation through a slate of statistics.

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

To write or read an instance of type *Colony* to or from a JSON file, call
```julia
toFile(colony, json_stream)
```
> which writes a colony of type *Colony* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
colony = fromFile(Colony, json_stream)
```
> reads a colony of type *Colony* from the JSON file attached to `json_stream`.
"""
struct Colony
    myparameters::AbstractParameters
    mydata::TheData
    # Fields pertaining to the creatures of a colony.
    probability_mutation::Float64
    probability_crossover::Float64
    probability_immigrant::Float64
    parameters_min::Vector{PF.PhysicalScalar}
    parameters_max::Vector{PF.PhysicalScalar}
    parameters_constrained::Vector{Tuple{Int,Int}}
    significant_figures::Int
    # Fields pertaining to the colony itself.
    population::Int
    generations::Int
    elite::Creature
    # Mutable fields pertaining to the colony itself.
    generation::PF.MInteger
    children::Vector{Creature}
    adults::Vector{Creature}

    # constructor

    function Colony(myparameters::AbstractParameters,
                    mydata::TheData, 
                    probability_mutation::Float64, 
                    probability_crossover::Float64, 
                    probability_immigrant::Float64,
                    parameters_alien::Vector{PF.PhysicalScalar},
                    parameters_min::Vector{PF.PhysicalScalar},
                    parameters_max::Vector{PF.PhysicalScalar}, 
                    parameters_constrained::Vector{Tuple{Int,Int}}, 
                    significant_figures::Int)
        
        # Bound the inputs.
        if significant_figures < 1
            significant_figures = 1
        end
        if significant_figures > 8
            significant_figures = 8
        end

        if probability_mutation < 1.0e-6
            probability_mutation = Float64(1.0e-6)
        end
        if probability_mutation > 0.999999
            probability_mutation = Float64(0.999999)
        end

        if probability_crossover < 1.0e-6
            probability_crossover = Float64(1.0e-6)
        end
        if probability_crossover > 0.999999
            probability_crossover = Float64(0.999999)
        end

        if probability_immigrant < 1.0e-6
            probability_immigrant = Float64(1.0e-6)
        end
        if probability_immigrant > 0.999999
            probability_immigrant = Float64(0.999999)
        end

        # Verify the dimensions.
        N = fieldcount(typeof(myparameters))

        if (length(parameters_alien) ≠ N) && (length(parameters_alien) ≠ 0)
            msg = "The supplied alien has the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end

        if (length(parameters_min) ≠ N) || (length(parameters_max) ≠ N)
            msg = "The supplied bounds have the wrong number of parameters."
            throw(DimensionMismatch(msg))
        end
        
        # Verify the units.     
        for n in 1:N
            if length(parameters_alien) == 0
                if parameters_min[n].units ≠ parameters_max[n].units
                    error("Parameters' physical units were not consistent.")
                end
            else
                if (parameters_min[n].units ≠ parameters_max[n].units ||
                    parameters_min[n].units ≠ parameters_alien[n].units)
                    error("Parameters' physical units were not consistent.")
                end
            end
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

        # Verify that the fields of object myparameters are PhysicalScalars.
        for n in 1:N
            symbol = fieldname(typeof(myparameters), n)
            if fieldtype(typeof(myparameters), symbol) ≠ PF.PhysicalScalar
                msg = "All parameter fields must be a PhysicalScalar.\n"
                msg = string(msg, "   Field ", n, " was a ")
                msg = string(msg, fieldtype(typeof(myparameters), symbol), ".")
                error(msg)
            end
        end

        # Create the model.
        mymodel = Model(myparameters, mydata)

        # Create the presupposed elite creature.
        elite = alien(parameters_alien, parameters_min, parameters_max,
                      parameters_constrained, significant_figures)
        set!(mymodel, parameters(elite))
        fitness = _evaluate(mymodel)
        PF.set!(elite.fitness, fitness)

        # Determine the population size.
        # Formula for estimating population size is from D. Goldberg (2002).
        alphabet = 2          # possible expressions: dominant and recessive
        schemata = significant_figures  # length of prominant gene sequences
        pop_size = (log(elite.DNA.genes) +
                    schemata * log(alphabet) * alphabet^schemata) 
        population = Int(ceil(pop_size))

        # Determine the required number of generations.
        # Formula for estimating convergence is from D. Goldberg (2002).
        combatants = population ÷ 50
        if combatants < 3
            combatants = 3
        end
        converge_at = sqrt(elite.DNA.genes) * log(population - combatants)
        generation  = PF.MInteger(1)
        generations = Int(ceil(converge_at))
       
        # Create the initial population of adults.
        adults = Vector{Creature}(undef, population)
        adults[1] = elite
        for i in 2:population
            adult = procreate(parameters_min, parameters_max,
                              parameters_constrained, significant_figures)
            thread_model = Model(myparameters, mydata)
            set!(thread_model, parameters(adult))
            thread_fitness = _evaluate(thread_model)
            PF.set!(adult.fitness, thread_fitness)
            adults[i] = adult
        end

        # Find the elite creature in the population.
        for i in 1:population
            adult = adults[i]
            if adult.fitness > elite.fitness
                elite = adult
            end
        end

        # Create the initial population of children.
        children = Vector{Creature}(undef, population)
        
        new(myparameters, mydata,
            probability_mutation, probability_crossover, probability_immigrant, 
            parameters_min, parameters_max, parameters_constrained,
            significant_figures, population, generations, elite, generation, 
            children, adults)::Colony
    end

    function Colony(myparameters::AbstractParameters,
                    mydata::TheData,
                    probability_mutation::Float64,
                    probability_crossover::Float64,
                    probability_immigrant::Float64,
                    parameters_min::Vector{PF.PhysicalScalar},
                    parameters_max::Vector{PF.PhysicalScalar},
                    parameters_constrained::Vector{Tuple{Int,Int}},
                    significant_figures::Int,
                    population::Int,
                    generations::Int,
                    elite::Creature,
                    generation::PF.MInteger,
                    children::Vector{Creature},
                    adults::Vector{Creature})

        new(myparameters, mydata,
            probability_mutation, probability_crossover, probability_immigrant, 
            parameters_min, parameters_max, parameters_constrained,
            significant_figures, population, generations, elite, generation, 
            children, adults)::Colony
    end
end # Colony

# internal methods

function _tournamentPlay(colony::Colony)::Creature
    if colony.probability_immigrant > rand()
        # Add an immigrant to the population.
        mostfit = procreate(colony.parameters_min, colony.parameters_max,
                            colony.parameters_constrained,
                            colony.significant_figures)
        mymodel = Model(colony.myparameters, colony.mydata)
        set!(mymodel, parameters(mostfit))
        fitness = _evaluate(mymodel)
        PF.set!(mostfit.fitness, fitness)
    else
        # Select a creature for mating.
        combatants = colony.population ÷ 50
        if combatants < 3
            combatants = 3
        end
        contestants = Vector{Creature}(undef, combatants)
        for i in 1:combatants
            combatant = rand(1:colony.population)
            contestants[i] = colony.adults[combatant]
        end
        mostfit = contestants[1]
        for combatant in 2:combatants
            contestant = contestants[combatant]
            if contestant.fitness > mostfit.fitness
                mostfit = contestant
            end
        end
    end
    return mostfit
end # _tournamentPlay

function _mate!(colony::Colony)
    # Elite creature from current generation lives into the next generation.
    colony.children[1] = copy(colony.elite)

    # Mate the population.
    for i in 2:colony.population
        parentA = _tournamentPlay(colony)
        parentB = _tournamentPlay(colony)
        self_conception = 0
        while parentB == parentA
            self_conception = self_conception + 1
            if self_conception == 25
                # This is for safety. It should not occur in practice.
                break
            end
            parentB = _tournamentPlay(colony)
        end
        child = conceive(parentA, parentB, 
                         colony.parameters_constrained,
                         colony.probability_mutation,
                         colony.probability_crossover)
        colony.children[i] = child
    end

    # Ensure there are no clones, i.e., that there are no identical twins.
    for i in 2:colony.population
        child = colony.children[i]
        for j in 1:i-1
            if child == colony.children[j]
                # Identical twins are not permitted within a generation.
                parentA = _tournamentPlay(colony)
                parentB = _tournamentPlay(colony)
                new_child = conceive(parentA, parentB,
                                     colony.parameters_constrained,
                                     colony.probability_mutation,
                                     colony.probability_crossover)
                colony.children[i] = new_child
                break
            end
        end
    end

    # Determine the fitness for each child born into the colony.
    for i in 2:colony.population
        child = colony.children[i]
        thread_model = Model(colony.myparameters, colony.mydata)
        set!(thread_model, parameters(child))
        fitness = _evaluate(thread_model)
        PF.set!(child.fitness, fitness)
        colony.children[i] = child
    end

    return nothing
end # _mate!

# Objective function

function _objfn(response_exp::PF.ArrayOfPhysicalScalars,
                response_mod::PF.ArrayOfPhysicalScalars)::Float64

    if response_exp.array.len == response_mod.array.len
        N = response_exp.array.len
    else
        error("Lengths of model and experimental response vectors must equal.")
    end
    if response_exp.units ≠ response_mod.units
        error("Physical units must equal for experiment and model responses.")
    end
    
    # Goal is to maximize an objective function => 1 / minimize the error.
    res_exp = response_exp.array.vec
    res_mod = response_mod.array.vec
    res_err = zeros(Float64, N)
    for n = 1:N
        res_err[n] = abs(res_exp[n] - res_mod[n])
    end
    
    objfn = (N / StatsBase.L1dist(res_exp, res_mod) +         # 1 / L₁ norm
             sqrt(N) / StatsBase.L2dist(res_exp, res_mod) +   # 1 / L₂ norm
             1 / StatsBase.Linfdist(res_exp, res_mod))        # 1 / L∞ norm

    return objfn
end # _objfn
    
function _evaluate(model::Model)::Float64

    # Populate the field model.data.response_mod.
    solve!(model)

    data_points = 0
    for exp in 1:model.data.experiments
        data_points += (model.data.conjugate_pairs[exp] *
                        model.data.data_points[exp])
    end

    fitness = 0.0
    for exp in 1:model.data.experiments
        ratio = model.data.data_points[exp] / data_points
        for pair in 1:model.data.conjugate_pairs[exp]
            ϕ = _objfn(model.data.response_exp[exp][pair],
                       model.data.response_mod[exp][pair])
            fitness += ratio * ϕ
        end
    end

    return fitness
end # _evaluate

function _corr_coef(model::Model)::Vector{Float64}
    data_points = 0
    for exp in 1:model.data.experiments
        data_points += (model.data.conjugate_pairs[exp] *
                        model.data.data_points[exp])
    end
    
    corr_coef = zeros(Float64, 4)
    for exp in 1:model.data.experiments
        ratio = model.data.data_points[exp] / data_points
        for pair in 1:model.data.conjugate_pairs[exp]
            res_exp = model.data.response_exp[exp][pair].array.vec
            res_mod = model.data.response_mod[exp][pair].array.vec
            R = Vector{Float64}(undef, 4)
            # Pearson's linear correlation coefficient r.
            R[1] = Statistics.cor(res_exp, res_mod)
            # Spearman's monotonic correlation coefficient ρ.
            R[2] = StatsBase.corspearman(res_exp, res_mod)
            # Kendall's monotonic correlation coefficient τ.
            R[3] = StatsBase.corkendall(res_exp, res_mod)
            # Chatterjee's nonlinear correlation coefficient ξ.
            R[4] = Xicor.xicor(res_exp, res_mod)
            for i in 1:4
                corr_coef[i] = corr_coef[i] + ratio * R[i]
            end
        end
    end
    
    return corr_coef
end # _corr_coef

function _stdDevElite(colony::Colony)::Vector{PF.PhysicalScalar}
    parameters_elite = parameters(colony.elite)
    len = length(parameters_elite)
    
    sum = zeros(Float64, len)
    for i in 1:colony.population
        adult = colony.adults[i]
        parameters_adult = parameters(adult)
        for j in 1:len
            sum[j] = sum[j] + (PF.get(parameters_adult[j]) -
                               PF.get(parameters_elite[j]))^2
        end
    end
    
    sde = Vector{PF.PhysicalScalar}(undef, len)
    for i in 1:len
        value  = sqrt(sum[i] / (colony.population - 1))
        units  = colony.elite.DNA.genotypes[i].parameter_min.units
        sde[i] = PF.PhysicalScalar(value, units)
    end
    
    return sde
end # _stdDevElite

# exported methods

function advance!(colony::Colony)
    PF.set!(colony.generation, PF.get(colony.generation)+1)
    _mate!(colony)

    # Children become the new adults.
    for i in 1:colony.population
        colony.adults[i] = copy(colony.children[i])
    end

    # Find the elite adult in this next generation.
    elite = colony.adults[1]
    for i in 2:colony.population
        adult = colony.adults[i]
        if adult.fitness > elite.fitness
            elite = adult
        end
    end
    PF.set!(colony.elite.fitness, PF.get(elite.fitness))
    for chromosome in 1:colony.elite.DNA.chromosomes
        colony.elite.DNA.genotypes[chromosome] = elite.DNA.genotypes[chromosome]
    end

    return nothing
end # advance!

function report(colony::Colony)::String
    g = PF.toString(colony.generation)
    p = string(colony.population)
    s = string("Fitness statistics for generation ", g)
    s = string(s, " with a population size of ", p)
    s = string(s, ":\n")
    fitness = PF.get(colony.elite.fitness)
    if fitness ≥ 0.0
        s = string(s, "   optimum fitness  ", PF.toString(fitness))
    else
        s = string(s, "   optimum fitness ", PF.toString(fitness))
    end
    s = string(s, ",\n")
    fitness_vec = Vector{Float64}(undef, colony.population)
    for i in 1:colony.population
        adult = colony.adults[i]
        fitness_vec[i] = get(adult.fitness)
    end
    stat_median = Statistics.median(fitness_vec)
    if stat_median ≥ 0.0
        s = string(s, "   median           ", PF.toString(stat_median))
    else
        s = string(s, "   median          ", PF.toString(stat_median))
    end
    s = string(s, ",\n")
    stat_mean = Statistics.mean(fitness_vec)
    if stat_mean ≥ 0.0
        s = string(s, "   arithmetic mean  ", PF.toString(stat_mean))
    else
        s = string(s, "   arithmetic mean ", PF.toString(stat_mean))
    end
    s = string(s, ",\n")
    stat_std = Statistics.std(fitness_vec)
    if stat_std ≥ 0.0
        s = string(s, "   std deviation    ", PF.toString(stat_std))
    else
        s = string(s, "   std deviation   ", PF.toString(stat_std))
    end
    s = string(s, ",\n")
    stat_skewness = StatsBase.skewness(fitness_vec)
    if stat_skewness ≥ 0.0
        s = string(s, "   skewness         ", PF.toString(stat_skewness))
    else
        s = string(s, "   skewness        ", PF.toString(stat_skewness))
    end
    s = string(s, ",\n")
    stat_kurtosis = StatsBase.kurtosis(fitness_vec)
    if stat_kurtosis ≥ 0.0
        s = string(s, "   excess kurtosis  ", PF.toString(stat_kurtosis))
    else
        s = string(s, "   excess kurtosis ", PF.toString(stat_kurtosis))
    end
    s = string(s, ".\n\n")
    if colony.generation == colony.generations
        mymodel = Model(colony.myparameters, colony.mydata)
        set!(mymodel, parameters(colony.elite))
        corr_coef = _corr_coef(mymodel)
        
        s = string(s, "To aid in assessing goodness of fit, ")
        s = string(s, "a suite of correlation coefficients")
        s = string(s, "\n")
        s = string(s, "between experimental response and model prediction ")
        s = string(s, "are provided")
        s = string(s, ".\n")
        s = string(s, "   Pearson's linear correlation coefficient")
        s = string(s, ":\n")
        s = string(s, "      r = ", PF.toString(corr_coef[1]), ",\n")
        s = string(s, "   Spearman's monotonic correlation coefficient")
        s = string(s, ":\n")
        s = string(s, "      ρ = ", PF.toString(corr_coef[2]), ",\n")
        s = string(s, "   Kendall's monotonic correlation coefficient")
        s = string(s, ":\n")
        s = string(s, "      τ = ", PF.toString(corr_coef[3]), ",\n")
        s = string(s, "   and Chatterjee's nonlinear correlation coefficient")
        s = string(s, ":\n")
        s = string(s, "      ξ = ", PF.toString(corr_coef[4]), ".\n\n")
    end
    s = string(s, "The genome from the most fit creature")
    s = string(s, ":\n")
    s = string(s, toBinaryString(colony.elite.DNA))
    s = string(s, "\n\n")
    s = string(s, "Lists for [parameter_min, parameter_best, parameter_max]")
    s = string(s, ":\n")
    parameters_best = parameters(colony.elite)
    
    # Get parameter names from the fields of object colony.myparameters.
    N = fieldcount(typeof(colony.myparameters))
    parameters_names = Vector{String}(undef, N)
    for n in 1:N
        symbol = fieldname(typeof(colony.myparameters), n)
        parameters_names[n] = chomp(String(symbol))
    end

    for n in 1:N
        if N < 10
            s = string(s, n, ": ", parameters_names[n], " ∈ ")
        else
            if n < 10
                s = string(s, n, ":  ", parameters_names[n], " ∈ ")
            else
                s = string(s, n, ": ", parameters_names[n], " ∈ ")
            end
        end
        s = string(s, "[", PF.toString(PF.get(colony.parameters_min[n])), ", ",
            PF.toString(PF.get(parameters_best[n])), ", ",
            PF.toString(PF.get(colony.parameters_max[n])), "] ",
            PF.toString(colony.parameters_min[n].units), "\n")
        if n == N
            s = string(s, "\n")
        end
    end
    s = string(s, "Values for best parameter ± error, where a RMSE ")
    s = string(s, "was computed wrt best values.")
    s = string(s, "\n")
    s = string(s, "Data were the parameters from all ", p, " adults living ")
    s = string(s, "in generation ", g)
    s = string(s, ".\n")
    err = _stdDevElite(colony)
    
    for n in 1:N
        if N < 10
            s = string(s, n, ": ", parameters_names[n], " = ")
        else
            if n < 10
                s = string(s, n, ":  ", parameters_names[n], " = ")
            else
                s = string(s, n, ": ", parameters_names[n], " = ")
            end
        end
        s = string(s, PF.toString(PF.get(parameters_best[n])), " ± ",
            PF.toString(err[n]), "\n")
    end
    s = string(s, "\n")
    return s
end # report

# Methods for storing and retrieving a Colony to and from a file.

StructTypes.StructType(::Type{Colony}) = StructTypes.Struct()

function toFile(colony::Colony, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, colony)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{Colony}, json_stream::IOStream)::Colony
    if isopen(json_stream)
        colony = JSON3.read(readline(json_stream), Colony)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return colony
end

