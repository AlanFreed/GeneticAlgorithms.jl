"""
A genome is a genetic container of chromosomes.  Here is where all the genetic
information exists that makes up a creature.

# Genome

```
struct Genome
    genes::Int64           # number of genes that comprise a genome
    chromosomes::Int64     # number of chromosomes that comprise a genome
    genotypes::Vector{Chromosome}
           # array of chromosomes, i.e., its genome
end
```
> where

Constructor

    g = Genome(parameters_min, parameters_max, parameters_units,
               parameters_constrained, significant_figures)

        parameters_min          array of most-negative parameter values θₘᵢₙ[i]
        parameters_max          array of most-positive parameter values θₘₐₓ[i]
        parameters_units        array of physical units for the parameters
        parameters_constrained  array of tuples (θL, θR) where θL[i] < θR[i]
        significant_figures     seek parameters with significant figure accuracy

Operators

    ==, ≠

Methods

    c = getindex(g, l)      return chromosome 'c' in genome 'g' at location 'l'
    setindex!(g, c, l)      assign chromosome 'c' to genome 'g' at location 'l'
    c = copy(g)             return a copy 'c' of genome 'g'
    s = toBinaryString(g)   return string 's' describing genome 'g'
    mutate!(g, pM)          random flip of gene expression with probability 'pM'
    crossover(A, B, pM, pX) crossover between chromosomes in genomes 'A' and 'B'
                            with probabilities of mutation 'pM' & crossover 'pX'
    θ = decode(g)           return phenotypes (the parameters 'θ') held by 'g'
    encode!(g, θ)           assign phenotypes (the parameters 'θ') to genome 'g'
"""
struct Genome
    genes::Int64
    chromosomes::Int64
    genotypes::Vector{Chromosome}

    # constructor

    function Genome(parameters_min::Vector{Float64}, parameters_max::Vector{Float64}, parameters_units::Vector{PhysicalFields.PhysicalUnits}, significant_figures::Int64)

        if length(parameters_min) == length(parameters_max)
            chromosomes = Int64(length(parameters_min))
        else
            msg = "Vectors parameters_min and parameters_max must have the same length."
            throw(DimensionMismatch, msg)
        end

        genotypes = Vector{Chromosome}(undef, chromosomes)
        for i in 1:chromosomes
            genotypes[i] = Chromosome(parameters_min[i], parameters_max[i], parameters_units[i], significant_figures)
        end

        genes = Int64(0)
        for i in 1:chromosomes
            chromosome = genotypes[i]
            genes = genes + chromosome.genes
        end

        new(genes, chromosomes, genotypes)::Genome
    end

    function Genome(parameters_min::Vector{Real}, parameters_max::Vector{Real}, parameters_units::Vector{PhysicalFields.PhysicalUnits}, significant_figures::Integer)

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

        return Genome(p_min, p_max, parameters_units, sigfig)
    end

    function Genome(genes::Int64, chromosomes::Int64, genotypes::Vector{Chromosome})

        new(genes, chromosomes, genotypes)::Genome
    end
end # Genome

# operators

function Base.:(==)(gL::Genome, gR::Genome)::Bool
    if gL.genes ≠ gR.genes
        return false
    end
    if gL.chromosomes ≠ gR.chromosomes
        return false
    end
    for i in 1:gL.chromosomes
        if gL.genotypes[i] ≠ gR.genotypes[i]
            return false
        end
    end
    return true
end # ==

function Base.:≠(gL::Genome, gR::Genome)::Bool
    if !(gL == gR)
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(getindex)(g::Genome, index::Int)::Chromosome
    if index ≥ 1 && index ≤ g.chromosomes
        chromosome = g.genotypes[index]
    else
        msg = string("Admissible chromosome indices are ∈ [1…", g.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return copy(chromosome)
end # getindex

function Base.:(setindex!)(g::Genome, chromosome::Chromosome, index::Int)
    if index ≥ 1 && index ≤ g.chromosomes
        old_chromosome = g.genotypes[index]
        if old_chromosome.parameter_units ≠ chromosome.parameter_units
            msg = "The chromosome to be assigned has the wrong physical units."
            error(msg)
        end
        g.genotypes[index] = copy(chromosome)
    else
        msg = string("Admissible chromosome indices are ∈ [1…", g.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return nothing
end # setindex!

function Base.:(copy)(g::Genome)::Genome
    genes       = copy(g.genes)
    chromosomes = copy(g.chromosomes)
    genotypes   = Vector{Chromosome}(undef, chromosomes)
    for i in 1:chromosomes
        genotypes[i] = copy(g.genotypes[i])
    end
    return Genome(genes, chromosomes, genotypes)
end # copy

function toBinaryString(g::Genome)::String
    s = ""
    for i in 1:g.chromosomes
        chromosome = g.genotypes[i]
        if g.chromosomes < 10
            s = string(s, i, ": ", toBinaryString(chromosome))
        else
            if i < 10
                s = string(s, i, ":  ", toBinaryString(chromosome))
            else
                s = string(s, i, ": ", toBinaryString(chromosome))
            end
        end
        if i < g.chromosomes
            s = string(s, "\n")
        end
    end
    return s
end # toBinaryString

function mutate!(g::Genome, probability_mutation::Float64)
    for i in 1:g.chromosomes
        mutate!(g.genotypes[i], probability_mutation)
    end
    return nothing
end # mutate!

function mutate!(g::Genome, probability_mutation::Real)
    prob_m = convert(Float64, probability_mutation)
    mutate!(g, prob_m)
    return nothing
end # mutate!

function crossover(parentA::Genome, parentB::Genome, probability_mutation::Float64, probability_crossover::Float64)::Genome

    if parentA.genes == parentB.genes
        if parentA.chromosomes == parentB.chromosomes
            child = copy(parentA)
            for i in 1:child.chromosomes
                child[i] = crossover(parentA[i], parentB[i], probability_mutation, probability_crossover)
            end
        else
            msg = "Crossover parents must have the same number of chromosomes."
            throw(DimensionMismatch(msg))
        end
    else
        msg = "Crossover parents must have the same number of genes."
        throw(DimensionMismatch(msg))
    end
    return child
end # crossover

function crossover(parentA::Genome, parentB::Genome, probability_mutation::Real, probability_crossover::Real)::Genome
    prob_m = convert(Float64, probability_mutation)
    prob_x = convert(Float64, probability_crossover)
    return crossover(parentA, parentB, prob_m, prob_x)
end # crossover

function decode(g::Genome)::Vector{PhysicalFields.PhysicalScalars}
    phenotypes = Vector{Float64}(undef, g.chromosomes)
    for i in 1:g.chromosomes
        chromosome = g.genotypes[i]
        phenotypes[i] = decode(chromosome)
    end
    return phenotypes
end # decode

function encode!(g::Genome, phenotypes::Vector{PhysicalFields.PhysicalScalars})
    if g.chromosomes ≠ length(phenotypes)
        msg = "Number of chromosomes must equal number of phenotypes."
        throw(DimensionMismatch, msg)
    end

    for i in 1:g.chromosomes
        chromosome = g.genotypes[i]
        phenotype  = phenotypes[i]
        encode!(chromosome, phenotype)
    end
end # encode!

# Methods for storing and retrieving a Genome to and from a file.

StructTypes.StructType(::Type{Genome}) = StructTypes.Struct()

function toFile(genome::Genome, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, genome)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

function fromFile(::Type{Genome}, json_stream::IOStream)::Genome
    if isopen(json_stream)
        genome = JSON3.read(readline(json_stream), Genome)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return genome
end
