"""
A chromosome is a genetic container of genes.  In this implementation of a
genetic algorithm, each chromosome represents a parameter (an unknown value in
some model to be parameterized).  Genetic processes (mutation and crossover)
adjust these chromosomes via evolution, and therefore, their associated
parametric values change over the generations.

Chromosomes are where genetics and optimization meet.

A chromosome has an interface of:

struct Chromosome
    parameter_min   # minimum value that a chromosome can represent
    parameter_max   # maximum value that a chromosome can represent
    genes           # number of genes that comprise a chromosome
    expressions     # number of gene expressions a chromosome can represent
    genotype        # genetic material, i.e., genes, comprising a chromosome
end

A parameter will be considered fixed if its parameter_min ≈ parameter_max.

Constructor

    c = Chromosome(parameter_min, parameter_max, significant_figures)

        parameter_min           least value a parameter can take on
        parameter_max           greatest value a parameter can take on
        significant_figures     seek parameter with significant figure accuracy

Operators

    ==, ≠

Methods

    g = getindex(c, i)      return gene 'g' from chromosome 'c' at index 'i'
    setindex!(c, g, i)      assign gene 'g' to chromosome 'c' at index 'i'
    d = copy(c)             return a copy 'd' of chromosome 'c'
    s = tostring(c)         return a string 's' describing chromosome 'c'
    θ = decode(c)           return a phenotype (the parameter 'θ') held by 'c'
    encode!(c, θ)           assign a phenotype 'θ' to chromosome 'c'
    mutate!(c, pM)          random flip of gene expressions at probability 'pM'
    crossover(A, B, pM, pX) crossover between chromosomes 'A' and 'B' at
                            probabilities of mutation 'pM' and crossover 'pX'
"""
struct Chromosome
    # Fields that bound a parameter.
    parameter_min::Real
    parameter_max::Real
    # Fields that describe a chromosome.
    genes::Integer
    expressions::Integer
    genotype::Vector{Gene}

    # constructor

    function Chromosome(parameter_min::Real, parameter_max::Real, significant_figures::Integer)

        if parameter_min ≈ parameter_max

            # Handle the case of a fixed parameter.
            parameter_min = parameter_max
            genes         = 0
            expressions   = 0
            genotype      = Vector{Gene}(undef, 0)

        else

            # Verify the input.
            if parameter_min > parameter_max
                msg = "Cannot create a chromosome unless parameter_max ≥ parameter_min."
                throw(ErrorException, msg)
            end

            # Number of decades that span [parameter_min, parameter_max].
            if parameter_min > 0.0
                logdecades = log10(parameter_max/parameter_min)
            elseif parameter_max < 0.0
                logdecades = log10(parameter_min/parameter_max)
            elseif parameter_min ≈ 0.0
                logdecades = log10(parameter_max)
            elseif parameter_max ≈ 0.0
                logdecades = log10(-parameter_min)
            else # parameter_min < 0.0 and parameter_max > 0.0
                logdecades = log10(-parameter_min*parameter_max)
            end
            if logdecades > 0.0
                decades = Int(ceil(logdecades))
            else
                decades = abs(Int(floor(logdecades)))
            end
            if decades < 1
                decades = 1
            end
            if decades > 9
                decades = 9
            end

            # Number of genes per decade of span in the parameter range.
            if significant_figures < 2
                genes_per_decade = 7
            elseif significant_figures == 2
                genes_per_decade = 10
            elseif significant_figures == 3
                genes_per_decade = 14
            elseif significant_figures == 4
                genes_per_decade = 17
            elseif significant_figures == 5
                genes_per_decade = 20
            elseif significant_figures == 6
                genes_per_decade = 24
            else
                genes_per_decade = 27
            end

            # Number of genes and the number of possible gene expressions.
            genes = genes_per_decade * decades
            if genes > 63
                genes = 63
            end
            if genes < 63
                expressions = 2^genes - 1
            else
                expressions = typemax(Int64)
            end

            # This constructor creates a random genotype
            genotype  = Vector{Gene}(undef, genes)
            if genes > 0
                # Use the following temporary variables.
                binary    = Vector{Bool}(undef, genes)
                gray      = Vector{Bool}(undef, genes)
                atinteger = rand(1:expressions)
                gene      = genes
                while atinteger > 0 && gene > 0
                    if atinteger % 2 == 0
                        binary[gene] = recessive
                    else
                        binary[gene] = dominant
                    end
                    atinteger = atinteger ÷ 2
                    gene      = gene - 1
                end
                # Remaining higher-order binary bits are zeros, i.e., recessive.
                for i in gene:-1:1
                    binary[i] = recessive
                end
                # Convert to a gray binary.
                gray[1] = binary[1]
                for i in 2:genes
                    gray[i] = binary[i-1] ⊻ binary[i]
                end
                # Assign this gray encoding to the genotype.
                for i in 1:genes
                    if gray[i] == dominant
                        genotype[i] = Gene(dominant)
                    else
                        genotype[i] = Gene(recessive)
                    end
                end
            end
        end

        new(parameter_min, parameter_max, genes, expressions, genotype)
    end
end # Chromosome

# operators

function Base.:(==)(cL::Chromosome, cR::Chromosome)::Bool
    if cL.genes ≠ cR.genes
        return false
    end
    if cL.genes == 0 && cL.parameter_min ≈ cR.parameter_min
        return true
    end
    for i in 1:cL.genes
        if cL.genotype[i] ≠ cR.genotype[i]
            return false
        end
    end
    return true
end # ==

function Base.:≠(cL::Chromosome, cR::Chromosome)::Bool
    if !(cL == cR)
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(getindex)(c::Chromosome, index::Integer)::Gene
    if c.genes == 0
        return nothing
    elseif index ≥ 1 && index ≤ c.genes
        gene = c.genotype[index]
    else
        msg = string("Admissible gene indices are ∈ [1…", c.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return copy(gene)
end # getindex

function Base.:(setindex!)(c::Chromosome, gene::Gene, index::Integer)
    if c.genes == 0
        return nothing
    elseif index ≥ 1 && index ≤ c.genes
        set!(c.genotype[index], get(gene))
    else
        msg = string("Admissible gene indices are ∈ [1…", c.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return nothing
end # setindex!

function Base.:(copy)(c::Chromosome)::Chromosome
    parameter_min = copy(c.parameter_min)
    parameter_max = copy(c.parameter_max)
    genes         = copy(c.genes)
    expressions   = copy(c.expressions)
    genotype      = Vector{Gene}(undef, genes)
    for i in 1:genes
        genotype[i] = copy(c.genotype[i])
    end
    return Chromosome(parameter_min, parameter_max, genes, expressions, genotype)
end # copy

function tostring(c::Chromosome)::String
    str = ""
    for i in 1:c.genes
        str = string(str, tostring(c.genotype[i]))
    end
    return str
end # tostring

function mutate!(c::Chromosome, probability_mutation::Real)
    for i in 1:c.genes
        mutate!(c.genotype[i], probability_mutation)
    end
    return nothing
end # mutate!

"""
Two chromosomes can reproduce creating an offspring, which is what procedure
crossover mimics.  It accepts two parent chromosomes along with probabilities of
occurrence for mutation and crossover (a splitting of the chromosome between a
random pair of neighboring genes).  The outcome is a child of the parents.
There is a small chance that the child will be a clone of one of its parents--a
chance that is dictated by the assigned probabilities.
"""
function crossover(parentA::Chromosome, parentB::Chromosome, probability_mutation::Real, probability_crossover::Real)::Chromosome

    if probability_crossover < 0.0 || probability_crossover ≥ 1.0
        msg = "A probability of crossover must belong to unit interval [0, 1)."
        throw(ErrorException, msg)
    end

    if parentA.genes == parentB.genes
        if rand() < 0.5
            child = copy(parentA)
            # Left side of chromosome splice belongs to parent A.
            if probability_crossover > rand()
                if child.genes > 5
                    xover = rand(3:child.genes-2)
                    for i in xover:child.genes
                        child[i] = parentB[i]
                    end
                elseif child.genes > 3
                    child[child.genes-1] = parentB[genes-1]
                    child[child.genes]   = parentB[genes]
                elseif child.genes > 0
                    child[child.genes] = parentB[genes]
                else
                    # do nothing
                end
            end
        else
            child = copy(parentB)
            # Left side of chromosome splice belongs to parent B.
            if probability_crossover > rand()
                if child.genes > 5
                    xover = rand(3:child.genes-2)
                    for i in xover:child.genes
                        child[i] = parentA[i]
                    end
                elseif child.genes > 3
                    child[child.genes-1] = parentA[genes-1]
                    child[child.genes]   = parentA[genes]
                elseif child.genes > 0
                    child[child.genes] = parentA[genes]
                else
                    # do nothing
                end
            end
        end
    else
        msg = "Crossover parents must have the same number of genes."
        throw(DimensionMismatch(msg))
    end

    mutate!(child, probability_mutation)

    return child
end # crossover

#=
The decode/encode map between a real phenotype and its haploid representation.

The algorithm for converting between binary and gray codes assumes the most
significant bit (MSB) is at the [Low] position of the code, while the least
significant bit (LSB) associates with the [High] position of the code, in other
words, e.g.,
    code = [1|0|1|1|0|0|1|0]
has a MSB of 1 and a LSB of 0.
=#

# Methods to decode a parameter.

function _gray2binary(gray::Vector{Bool})::Vector{Bool}
    bits   = length(gray)
    binary = Vector{Bool}(undef, bits)
    binary[1] = gray[1]
    for i in 2:bits
        binary[i] = binary[i-1] ⊻ gray[i]
    end
    return binary
end # _gray2binary

function _binary2integer(binary::Vector{Bool})::Integer
    bits    = length(binary)
    integer = 0
    power   = 1
    for i in bits:-1:1
        if binary[i] == dominant
            integer = integer + power
        end
        power = 2power
    end
    return integer
end # _binary2integer

function _integer2phenotype(c::Chromosome, integer::Integer)::Real
    phenotype = (c.parameter_min + (integer / c.expressions)
        * (c.parameter_max - c.parameter_min))
    if phenotype < c.parameter_min
        phenotype = c.parameter_min
    end
    if phenotype > c.parameter_max
        phenotype = c.parameter_max
    end
    return phenotype
end # _integer2phenotype

function decode(c::Chromosome)::Real
    if c.genes == 0
        phenotype = copy(c.parameter_min)
    else
        gray = Vector{Bool}(undef, c.genes)
        for i in 1:c.genes
            if isdominant(c.genotype[i])
                gray[i] = dominant
            else
                gray[i] = recessive
            end
        end
        binary    = _gray2binary(gray)
        integer   = _binary2integer(binary)
        phenotype = _integer2phenotype(c, integer)
    end
    return phenotype
end # decode

# Methods to encode a parameter.

function _phenotype2integer(c::Chromosome, phenotype::Real)::Integer
    fraction = (phenotype - c.parameter_min) / (c.parameter_max - c.parameter_min)
    integer  = Int(round(fraction * c.expressions))
    if integer < 0
        integer = 0
    end
    if integer > typemax(Int64)
        integer = typemax(Int64)
    end
    return integer
end # _phenotype2integer

function _integer2binary(c::Chromosome, integer::Integer)::Vector{Bool}
    binary    = Vector{Bool}(undef, c.genes)
    atinteger = integer
    gene      = c.genes
    while atinteger > 0
        if atinteger % 2 == 0
            binary[gene] = recessive
        else
            binary[gene] = dominant
        end
        atinteger = atinteger ÷ 2
        gene      = gene - 1
    end
    # The remaining higher-order binary bits are zeros, i.e., recessive genes.
    for i in gene:-1:1
        binary[i] = recessive
    end
    return binary
end # _integer2binary

function _binary2gray(binary::Vector{Bool})::Vector{Bool}
    bits = length(binary)
    gray = Vector{Bool}(undef, bits)
    gray[1] = binary[1]
    for i in 2:bits
        gray[i] = binary[i-1] ⊻ binary[i]
    end
    return gray
end # _binary2gray

function encode!(c::Chromosome, phenotype::Real)
    if c.genes == 0
        # do nothing
    elseif phenotype ≥ c.parameter_min && phenotype ≤ c.parameter_max
        integer = _phenotype2integer(c, phenotype)
        binary  = _integer2binary(c, integer)
        gray    = _binary2gray(binary)
        for i in 1:c.genes
            if gray[i] == dominant
                c.genotype[i] = Gene(dominant)
            else
                c.genotype[i] = Gene(recessive)
            end
        end
    else
        msg = "The phenotype must lie within [c.parameter_min, c.parameter_max]."
        throw(ErrorException, msg)
    end
    return nothing
end # encode
