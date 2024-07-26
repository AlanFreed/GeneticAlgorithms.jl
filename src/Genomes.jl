"""
A genome is a genetic container of chromosomes.  Here is where all the genetic
information exists that makes up a creature.

A genome has an interface of:

struct Genome
    genes           # number of genes that comprise a genome
    chromosomes     # number of chromosomes that comprise a genome
    genotypes       # array of chromosomes, i.e., its genome
end

Constructor

    g = Genome(parameters_min, parameters_max, parameters_constrained,
               significant_figures)

        parameters_min          array of most-negative parameter values θₘᵢₙ[i]
        parameters_max          array of most-positive parameter values θₘₐₓ[i]
        parameters_constrained  array of tuples (θL, θR) where θL[i] < θR[i]
        significant_figures     seek parameters with significant figure accuracy

Operators

    ==, ≠

Methods

    c = getindex(g, l)      return chromosome 'c' in genome 'g' at location 'l'
    setindex!(g, c, l)      assign chromosome 'c' to genome 'g' at location 'l'
    c = copy(g)             return a copy 'c' of genome 'g'
    s = tostring(g)         return string 's' describing genome 'g'
    mutate!(g, pM)          random flip of gene expression with probability 'pM'
    crossover(A, B, pM, pX) crossover between chromosomes in genomes 'A' and 'B'
                            with probabilities of mutation 'pM' & crossover 'pX'
    θ = decode(g)           return phenotypes (the parameters 'θ') held by 'g'
    encode!(g, θ)           assign phenotypes (the parameters 'θ') to genome 'g'
"""
struct Genome
    genes::Integer
    chromosomes::Integer
    genotypes::Vector{Chromosome}

    # constructor

    function Genome(parameters_min::Vector{Real}, parameters_max::Vector{Real}, significant_figures::Integer)

        if length(parameters_min) == length(parameters_max)
            chromosomes = length(parameters_min)
        else
            msg = "Vectors parameters_min and parameters_max must have the same length."
            throw(DimensionMismatch, msg)
        end

        genotypes = Vector{Chromosome}(undef, chromosomes)
        for i in 1:chromosomes
            genotypes[i] = Chromosome(parameters_min[i], parameters_max[i], significant_figures)
        end

        genes = 0
        for i in 1:chromosomes
            chromosome = genotypes[i]
            genes = genes + chromosome.genes
        end

        new(genes, chromosomes, genotypes)
    end

    function Genome(genes::Integer, chromosomes::Integer, genotypes::Vector{Chromosome})

        new(genes, chromosomes, genotypes)
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

function Base.:(getindex)(g::Genome, index::Integer)::Chromosome
    if index ≥ 1 && index ≤ g.chromosomes
        chromosome = g.genotypes[index]
    else
        msg = string("Admissible chromosome indices are ∈ [1…", g.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return copy(chromosome)
end # getindex

function Base.:(setindex!)(g::Genome, chromosome::Chromosome, index::Integer)
    if index ≥ 1 && index ≤ g.chromosomes
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

function tostring(g::Genome)::String
    str = ""
    for i in 1:g.chromosomes
        if g.chromosomes < 10
            str = string(str, i, ": ", tostring(g.genotypes[i]))
        else
            if i < 10
                str = string(str, i, ":  ", tostring(g.genotypes[i]))
            else
                str = string(str, i, ": ", tostring(g.genotypes[i]))
            end
        end
        if i < g.chromosomes
            str = string(str, "\n")
        end
    end
    return str
end # tostring

function mutate!(g::Genome, probability_mutation::Real)
    for i in 1:g.chromosomes
        mutate!(g.genotypes[i], probability_mutation)
    end
end # mutate!

function crossover(parentA::Genome, parentB::Genome, probability_mutation::Real, probability_crossover::Real)::Genome

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

function decode(g::Genome)::Vector{Real}
    phenotypes = Vector{Real}(undef, g.chromosomes)
    for i in 1:g.chromosomes
        chromosome = g.genotypes[i]
        phenotypes[i] = decode(chromosome)
    end
    return phenotypes
end # decode

function encode!(g::Genome, phenotypes::Vector{Real})
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
