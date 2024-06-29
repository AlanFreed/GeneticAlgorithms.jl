"""
A genome is a genetic container of chromosomes.  Here is where all the genetic
information exists that makes a creature unique from other creatures.

A genome has an interface of:

struct Genome
    genes           # number of genes that comprise a genome
    chromosomes     # number of chromosomes that comprise a genome
    genotypes       # an array of chromosomes, i.e., its genome
end

Constructors

    g = Genome(minParameters, maxParameters, significantFigures)

        minParameters       array of most-negative values parameters can have
        maxParameters       array of most-positive values parameters can have
        significantFigures  seek parameters with significant figure accuracy
or
    g = Genome(genes, chromosomes, genotypes)

Operators

    ==, ≠

Methods

    c = getindex(g, l)      return chromosome 'c' in genome 'g' at location 'l'
    setindex!(g, c, l)      assign chromosome 'c' to genome 'g' at location 'l'
    c = copy(g)             return a copy 'c' of genome 'g'
    c = deepcopy(g)         return a deep copy 'c' of genome 'g'
    s = toString(g)         return string 's' describing genome 'g'
    mutate(g, probability)  random flip of gene expression with 'probability'
    crossover(A, B, pM, pX) crossover between chromosomes in genomes 'A' and 'B'
                            with probabilities of mutation 'pM' & crossover 'pX'
    p = decode(g)           return phenotypes (the parameters 'p') held by 'g'
    encode!(g, p)           assign phenotypes 'p' to genome 'g'
"""
struct Genome
    genes::Int64
    chromosomes::Int64
    genotypes::Vector{Chromosome}

    # constructors

    function Genome(minParameters::Vector{Float64}, maxParameters::Vector{Float64}, significantFigures::Int64)

        if length(minParameters) == length(maxParameters)
            chromosomes = length(minParameters)
        else
            msg = "Vectors minParameters and maxParameters must have the same length."
            throw(DimensionMismatch, msg)
        end

        genotypes = Vector{Chromosome}(undef, chromosomes)
        for i in 1:chromosomes
            genotypes[i] = Chromosome(minParameters[i], maxParameters[i], significantFigures)
        end

        genes = 0
        for i in 1:chromosomes
            genes = genes + genotypes[i].genes
        end

        new(genes, chromosomes, genotypes)
    end

    function Genome(genes::Int64, chromosomes::Int64, genotypes::Vector{Chromosome})
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
    if (index ≥ 1) && (index ≤ g.chromosomes)
        chromosome = g.genotypes[index]
    else
        msg = string("Admissible chromosome indices are ∈ [1…", g.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return deepcopy(chromosome)
end # getindex

function Base.:(setindex!)(g::Genome, chromosome::Chromosome, index::Integer)
    if (index ≥ 1) && (index ≤ g.chromosomes)
        g.genotypes[index] = deepcopy(chromosome)
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

function Base.:(deepcopy)(g::Genome)::Genome
    genes       = deepcopy(g.genes)
    chromosomes = deepcopy(g.chromosomes)
    genotypes   = Vector{Chromosome}(undef, chromosomes)
    for i in 1:chromosomes
        genotypes[i] = deepcopy(g.genotypes[i])
    end
    return Genome(genes, chromosomes, genotypes)
end # deepcopy

function toString(g::Genome)::String
    str = ""
    for i in 1:g.chromosomes
        if g.chromosomes < 10
            str = string(str, i, ": ", toString(g.genotypes[i]))
        else
            if i < 10
                str = string(str, i, ":  ", toString(g.genotypes[i]))
            else
                str = string(str, i, ": ", toString(g.genotypes[i]))
            end
        end
        if i < g.chromosomes
            str = string(str, "\n")
        end
    end
    return str
end # toString

function mutate!(g::Genome, probabilityOfMutation::Float64)
    for i in 1:g.chromosomes
        mutate!(g.genotypes[i], probabilityOfMutation)
    end
end # mutate!

function crossover(parentA::Genome, parentB::Genome, probabilityOfMutation::Float64, probabilityOfCrossover::Float64)::Genome

    if parentA.genes == parentB.genes
        if parentA.chromosomes == parentB.chromosomes
            child = deepcopy(parentA)
            for i in 1:child.chromosomes
                child[i] = crossover(parentA[i], parentB[i], probabilityOfMutation, probabilityOfCrossover)
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

function decode(g::Genome)::Vector{Float64}
    phenotypes = Vector{Float64}(undef, g.chromosomes)
    for i in 1:g.chromosomes
        chromosome = g.genotypes[i]
        phenotypes[i] = decode(chromosome)
    end
    return phenotypes
end # decode

function encode!(g::Genome, phenotypes::Vector{Float64})
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
