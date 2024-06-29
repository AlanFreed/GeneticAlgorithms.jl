"""
A chromosome is a genetic container of genes.  In this implementation of a
genetic algorithm, each chromosome represents a parameter (an unknown value in
some model to be parameterized).  Genetic processes (mutation and crossover)
adjust these chromosomes via evolution, and therefore, their associated
parametric values change over the generations.

Chromosomes are where genetics and optimization meet.

A chromosome has an interface of:

struct Chromosome
    minParameter    # minimum value that a chromosome can represent
    maxParameter    # maximum value that a chromosome can represent
    genes           # number of genes that comprise a chromosome
    expressions     # number of gene expressions a chromosome can represent
    genotype        # genetic material, i.e., genes, comprising a chromosome
end

A parameter will be considered fixed if its minParameter ≈ maxParameter.

Constructors

    c = Chromosome(minParameter, maxParameter, significantFigures)

        minParameter        most negative value a parameter can take on
        maxParameter        most positive value a parameter can take on
        significantFigures  seek parameter with significant figure accuracy
or
    c = Chromosome(minParameter, maxParameter, genes, expressions, genotype)

Operators

    ==, ≠

Methods

    g = getindex(c, l)      return gene 'g' from chromosome 'c' at index 'l'
    setindex!(c, g, l)      assign gene 'g' to chromosome 'c' at location 'l'
    d = copy(c)             return a copy 'd' of chromosome 'c'
    d = deepcopy(c)         return a deep copy 'd' of chromosome 'c'
    s = toString(c)         return a string 's' describing chromosome 'c'
    p = decode(c)           return a phenotype (the parameter 'p') held by 'c'
    encode!(c, p)           assigns a phenotype 'p' to chromosome 'c'
    mutate(c, probability)  random flip of gene expressions at a 'probability'
    crossover(A, B, pM, pX) crossover between chromosomes 'A' and 'B' at
                            probabilities of mutation 'pM' and crossover 'pX'
"""
struct Chromosome
    # Fields that bound a parameter.
    minParameter::Float64
    maxParameter::Float64
    # Fields that describe a chromosome.
    genes::Int64
    expressions::Int64
    genotype::Vector{Gene}

    # constructor

    function Chromosome(minParameter::Float64, maxParameter::Float64, significantFigures::Int64)

        if minParameter ≈ maxParameter

            # Handle the case of a fixed parameter.
            minParameter = maxParameter
            genes        = 0
            expressions  = 0 
            genotype     = Vector{Gene}(undef, 0)

        else

            # Verify the input.
            if minParameter > maxParameter
                msg = "Cannot create chromosome unless maxParameter > minParameter."
                throw(ErrorException, msg)
            end

            # Number of decades that span [minParameter, maxParameter].
            if minParameter > 0.0
                logDecades = log10(maxParameter/minParameter)
            elseif maxParameter < 0.0
                logDecades = log10(minParameter/maxParameter)
            elseif minParameter == 0.0
                logDecades = log10(maxParameter)
            elseif maxParameter == 0.0
                logDecades = log10(-minParameter)
            else # minParameter < 0.0 and maxParameter > 0.0
                logDecades = log10(-minParameter*maxParameter)
            end
            if logDecades > 0.0
                decades = Int64(ceil(logDecades))
            else
                decades = abs(Int64(floor(logDecades)))
            end
            if decades < 1
                decades = 1
            end
            if decades > 9
                decades = 9
            end

            # Number of genes per decade of span in the parameter range.
            if significantFigures < 2
                genesPerDecade = 7
            elseif significantFigures == 2
                genesPerDecade = 10
            elseif significantFigures == 3
                genesPerDecade = 14
            elseif significantFigures == 4
                genesPerDecade = 17
            elseif significantFigures == 5
                genesPerDecade = 20
            elseif significantFigures == 6
                genesPerDecade = 24
            else
                genesPerDecade = 27
            end

            # Number of genes and the number of possible gene expressions.
            genes = genesPerDecade * decades
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
            # using the following temporary variables.
            binary    = Vector{Bool}(undef, genes)
            gray      = Vector{Bool}(undef, genes)
            atInteger = Int64(rand(1:expressions))
            gene      = genes
            while (atInteger > 0) && (gene > 0)
                if atInteger % 2 == 0
                    binary[gene] = recessive
                else
                    binary[gene] = dominant
                end
                atInteger = atInteger ÷ 2
                gene      = gene - 1
            end
            # The remaining higher-order binary bits are zeros, i.e., recessive.
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

        new(minParameter, maxParameter, genes, expressions, genotype)
    end

    function Chromosome(minParameter::Float64, maxParameter::Float64, genes::Int64, expressions::Int64, genotype::Vector{Gene})

        new(minParameter, maxParameter, genes, expressions, genotype)
    end
end # Chromosome

# operators

function Base.:(==)(cL::Chromosome, cR::Chromosome)::Bool
    if cL.genes ≠ cR.genes
        return false
    end
    if (cL.genes == 0) && (cL.minParameter ≈ cR.minParameter)
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
    elseif (index ≥ 1) && (index ≤ c.genes)
        gene = c.genotype[index]
    else
        msg = string("Admissible gene indices are ∈ [1…", c.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return deepcopy(gene)
end # getindex

function Base.:(setindex!)(c::Chromosome, gene::Gene, index::Integer)
    if c.genes == 0
        return nothing
    elseif (index ≥ 1) && (index ≤ c.genes)
        set!(c.genotype[index], get(gene))
    else
        msg = string("Admissible gene indices are ∈ [1…", c.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return nothing
end # setindex!

function Base.:(copy)(c::Chromosome)::Chromosome
    minParameter = copy(c.minParameter)
    maxParameter = copy(c.maxParameter)
    genes        = copy(c.genes)
    expressions  = copy(c.expressions)
    genotype     = Vector{Gene}(undef, genes)
    for i in 1:genes
        genotype[i] = copy(c.genotype[i])
    end
    return Chromosome(minParameter, maxParameter, genes, expressions, genotype)
end # copy

function Base.:(deepcopy)(c::Chromosome)::Chromosome
    minParameter = deepcopy(c.minParameter)
    maxParameter = deepcopy(c.maxParameter)
    genes        = deepcopy(c.genes)
    expressions  = deepcopy(c.expressions)
    genotype     = Vector{Gene}(undef, genes)
    for i in 1:genes
        genotype[i] = deepcopy(c.genotype[i])
    end
    return Chromosome(minParameter, maxParameter, genes, expressions, genotype)
end # deepcopy

function toString(c::Chromosome)::String
    str = ""
    for i in 1:c.genes
        str = string(str, toString(c.genotype[i]))
    end
    return str
end # toString

function mutate!(c::Chromosome, probabilityOfMutation::Float64)
    for i in 1:c.genes
        mutate!(c.genotype[i], probabilityOfMutation)
    end
end # mutate!

"""
Two chromosomes can reproduce creating an offspring, which is what procedure
crossover mimics.  It accepts two parent chromosomes along with probabilities of
occurrence for mutation and crossover (a splitting of the chromosome between a
random pair of neighboring genes).  The outcome is a child of the parents.
There is a small chance that the child will be a clone of one of its parents--a
chance that is dictated by the assigned probabilities.
"""
function crossover(parentA::Chromosome, parentB::Chromosome, probabilityOfMutation::Float64, probabilityOfCrossover::Float64)::Chromosome

    if (probabilityOfCrossover < 0.0) || (probabilityOfCrossover ≥ 1.0)
        msg = "A probability of crossover must belong to unit interval [0, 1)."
        throw(ErrorException, msg)
    end

    if parentA.genes == parentB.genes
        if rand() < 0.5
            child = deepcopy(parentA)
            # Left side of chromosome splice belongs to parent A.
            if probabilityOfCrossover > rand()
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
            child = deepcopy(parentB)
            # Left side of chromosome splice belongs to parent B.
            if probabilityOfCrossover > rand()
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

    mutate!(child, probabilityOfMutation)

    return child
end # crossover

#=
The decode/encode maps between a real phenotype and its Haploid representation.

The algorithm for converting between binary and gray codes assumes the most
significant bit (MSB) is at the [Low] position of the code, while the least
significant bit (LSB) associates with the [High] position of the code, in other
words, e.g.,
    code = [1|0|1|1|0|0|1|0]
has a MSB of 1 and a LSB of 0.
=#

# Methods to decode a parameter.

function _greyToBinary(grey::Vector{Bool})::Vector{Bool}
    bits   = length(grey)
    binary = Vector{Bool}(undef, bits)
    binary[1] = grey[1]
    for i in 2:bits
        binary[i] = binary[i-1] ⊻ grey[i]
    end
    return binary
end # _greyToBinary

function _binaryToInteger(binary::Vector{Bool})::Int64
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
end # _binaryToInteger

function _integerToPhenotype(c::Chromosome, integer::Int64)::Float64
    phenotype = (c.minParameter + (Float64(integer) / Float64(c.expressions))
        * (c.maxParameter - c.minParameter))
    if phenotype < c.minParameter
        phenotype = c.minParameter
    end
    if phenotype > c.maxParameter
        phenotype = c.maxParameter
    end
    return phenotype
end # _integerToPhenotype

function decode(c::Chromosome)::Float64
    if c.genes == 0
        phenotype = deepcopy(c.minParameter)
    else
        grey = Vector{Bool}(undef, c.genes)
        for i in 1:c.genes
            if isDominant(c.genotype[i])
                grey[i] = dominant
            else
                grey[i] = recessive
            end
        end
        binary    = _greyToBinary(grey)
        integer   = _binaryToInteger(binary)
        phenotype = _integerToPhenotype(c, integer)
    end
    return phenotype
end # decode

# Methods to encode a parameter.

function _phenotypeToInteger(c::Chromosome, phenotype::Float64)::Int64
    fraction = (phenotype - c.minParameter) / (c.maxParameter - c.minParameter)
    integer  = Int64(round(fraction * Float64(c.expressions)))
    if integer < 0
        integer = 0
    end
    if integer > typemax(Int64)
        integer = typemax(Int64)
    end
    return integer
end # _phenotypeToInteger

function _integerToBinary(c::Chromosome, integer::Int64)::Vector{Bool}
    binary = Vector{Bool}(undef, c.genes)
    atInt  = integer
    gene   = c.genes
    while atInt > 0
        if atInt % 2 == 0
            binary[gene] = recessive
        else
            binary[gene] = dominant
        end
        atInt = atInt ÷ 2
        gene  = gene - 1
    end
    # The remaining higher-order binary bits are zeros, i.e., recessive.
    for i in gene:-1:1
        binary[i] = recessive
    end
    return binary
end # _integerToBinary

function _binaryToGrey(binary::Vector{Bool})::Vector{Bool}
    bits = length(binary)
    grey = Vector{Bool}(undef, bits)
    grey[1] = binary[1]
    for i in 2:bits
        grey[i] = binary[i-1] ⊻ binary[i]
    end
    return grey
end # _binaryToGrey

function encode!(c::Chromosome, phenotype::Float64)
    if c.genes == 0
        # do nothing
    elseif (phenotype ≥ c.minParameter) && (phenotype ≤ c.maxParameter)
        integer = _phenotypeToInteger(c, phenotype)
        binary  = _integerToBinary(c, integer)
        grey    = _binaryToGrey(binary)
        for i in 1:c.genes
            if grey[i] == dominant
                c.genotype[i] = Gene(dominant)
            else
                c.genotype[i] = Gene(recessive)
            end
        end
    else
        msg = "The phenotype must lie within [c.minParameter, c.maxParameter]."
        throw(ErrorException, msg)
    end
    return nothing
end # encode
