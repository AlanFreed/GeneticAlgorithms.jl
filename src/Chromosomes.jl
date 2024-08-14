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
    parameter_units # physical units of the parameter
    genes           # number of genes that comprise a chromosome
    expressions     # number of gene expressions a chromosome can represent
    genotype        # genetic material, i.e., genes, comprising a chromosome
end

A parameter will be considered fixed if its parameter_min ≈ parameter_max.

Constructor

    c = Chromosome(parameter_min, parameter_max, parameter_units,
                   significant_figures)

        parameter_min           least value a parameter can take on
        parameter_max           greatest value a parameter can take on
        parameter_units         physical units of the parameter
        significant_figures     seek parameter with significant figure accuracy

Operators

    ==, ≠

Methods

    g = getindex(c, i)      return gene 'g' from chromosome 'c' at index 'i'
    setindex!(c, g, i)      assign gene 'g' to chromosome 'c' at index 'i'
    d = copy(c)             return a copy 'd' of chromosome 'c'
    s = toBinaryString(c)   return a string 's' describing chromosome 'c'
    θ = decode(c)           return a phenotype (the parameter 'θ') held by 'c'
    encode!(c, θ)           assign a phenotype 'θ' to chromosome 'c'
    mutate!(c, pM)          random flip of gene expressions at probability 'pM'
    crossover(A, B, pM, pX) crossover between chromosomes 'A' and 'B' at
                            probabilities of mutation 'pM' and crossover 'pX'
"""
struct Chromosome
    # Fields that bound a parameter.
    parameter_min::Float64
    parameter_max::Float64
    parameter_units::PhysicalFields.PhysicalUnits
    # Fields that describe a chromosome.
    genes::Int64
    expressions::Int64
    genotype::Vector{Gene}

    # constructor

    function Chromosome(parameter_min::Float64, parameter_max::Float64, parameter_units::PhysicalFields.PhysicalUnits, significant_figures::Int64)

        if parameter_min ≈ parameter_max

            # Handle the case of a fixed parameter.
            parameter_min = parameter_max
            genes         = Int64(0)
            expressions   = Int64(0)
            genotype      = Vector{Gene}(undef, 0)

        else

            # Verify the input.
            if parameter_min > parameter_max
                msg = "Cannot create a chromosome unless parameter_max ≥ parameter_min."
                error(msg)
            end

            # Number of decades that span [parameter_min, parameter_max].
            if parameter_min > 0.0
                log_decades = log10(parameter_max/parameter_min)
            elseif parameter_max < 0.0
                log_decades = log10(parameter_min/parameter_max)
            elseif parameter_min ≈ 0.0
                log_decades = log10(parameter_max)
            elseif parameter_max ≈ 0.0
                log_decades = log10(-parameter_min)
            else # parameter_min < 0.0 and parameter_max > 0.0
                log_decades = log10(-parameter_min*parameter_max)
            end
            if log_decades > 0.0
                decades = Int(ceil(log_decades))
            else
                decades = Int(abs(floor(log_decades)))
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
            genes = Int64(genes_per_decade * decades)
            if genes > 63
                genes = Int64(63)
            end
            if genes < 63
                expressions = Int64(2^genes - 1)
            else
                expressions = typemax(Int64)
            end

            # This constructor creates a random genotype
            genotype = Vector{Gene}(undef, genes)
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

        new(parameter_min, parameter_max, parameter_units, genes, expressions, genotype)::Chromosome
    end

    function Chromosome(parameter_min::Real, parameter_max::Real, parameter_units::PhysicalFields.PhysicalUnits, significant_figures::Integer)

        p_min  = convert(Float64, parameter_min)
        p_max  = convert(Float64, parameter_max)
        sigfig = convert(Int64, significant_figures)

        return Chromosome(p_min, p_max, parameter_units, sigfig)
    end

    function Chromosome(parameter_min::Float64, parameter_max::Float64, parameter_units::PhysicalFields.PhysicalUnits, genes::Int64, expressions::Int64, genotype::Vector{Gene})

        new(parameter_min, parameter_max, parameter_units, genes, expressions, genotype)::Chromosome
    end
end # Chromosome

# operators

function Base.:(==)(cL::Chromosome, cR::Chromosome)::Bool
    if cL.parameter_units ≠ cR.parameter_units || cL.genes ≠ cR.genes
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

function Base.:(getindex)(c::Chromosome, index::Int)::Gene
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

function Base.:(setindex!)(c::Chromosome, gene::Gene, index::Int)
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
    parameter_min   = copy(c.parameter_min)
    parameter_max   = copy(c.parameter_max)
    parameter_units = PhysicalFields.copy(c.parameter_units)
    genes           = copy(c.genes)
    expressions     = copy(c.expressions)
    genotype        = Vector{Gene}(undef, genes)
    for i in 1:genes
        genotype[i] = copy(c.genotype[i])
    end
    return Chromosome(parameter_min, parameter_max, parameter_units, genes, expressions, genotype)
end # copy

function toBinaryString(c::Chromosome)::String
    s = ""
    for i in 1:c.genes
        gene = c.genotype[i]
        s = string(s, toBinaryString(gene))
    end
    return s
end # toBinaryString

function mutate!(c::Chromosome, probability_mutation::Float64)
    for i in 1:c.genes
        gene = c.genotype[i]
        mutate!(gene, probability_mutation)
    end
    return nothing
end # mutate!

function mutate!(c::Chromosome, probability_mutation::Real)
    probability = convert(Float64, probability_mutation)
    for i in 1:c.genes
        gene = c.genotype[i]
        mutate!(gene, probability)
    end
    return nothing
end # mutate!

"""
Two chromosomes can reproduce creating an offspring, which is what procedure
crossover mimics.  It accepts two parent chromosomes along with probabilities of
occurrence for mutation (a gene swapping its expression) and crossover (a
splitting of the chromosome between a random pair of neighboring genes).  The
outcome is a child of the parents.  There is a small chance that the child will
be a clone of one of its parents--a chance that is dictated by the assigned
probabilities.
"""
function crossover(parentA::Chromosome, parentB::Chromosome, probability_mutation::Float64, probability_crossover::Float64)::Chromosome

    if probability_crossover < 0.0 || probability_crossover ≥ 1.0
        msg = "A probability of crossover must belong to unit interval [0, 1)."
        error(msg)
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

function crossover(parentA::Chromosome, parentB::Chromosome, probability_mutation::Real, probability_crossover::Real)::Chromosome
    probability_m = convert(Float64, probability_mutation)
    probability_c = convert(Float64, probability_crossover)
    child = crossover(parentA, parentB, probability_m, probabilty_c)
    return child
end # crossover

#=
The decode/encode map between a real phenotype and its haploid genotype
representation.

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

function _binary2integer(binary::Vector{Bool})::Int64
    bits    = length(binary)
    integer = 0
    power   = 1
    for i in bits:-1:1
        if binary[i] == dominant
            integer = integer + power
        end
        power = 2power
    end
    return convert(Int64, integer)
end # _binary2integer

function _integer2phenotype(c::Chromosome, integer::Int64)::Float64
    phenotype = (c.parameter_min + (integer / c.expressions)
        * (c.parameter_max - c.parameter_min))
    if phenotype < c.parameter_min
        phenotype = c.parameter_min
    end
    if phenotype > c.parameter_max
        phenotype = c.parameter_max
    end
    return convert(Float64, phenotype)
end # _integer2phenotype

function decode(c::Chromosome)::PhysicalFields.PhysicalScalar
    if c.genes == 0
        phenotype = copy(c.parameter_min)
    else
        gray = Vector{Bool}(undef, c.genes)
        for i in 1:c.genes
            if isDominant(c.genotype[i])
                gray[i] = dominant
            else
                gray[i] = recessive
            end
        end
        binary    = _gray2binary(gray)
        integer   = _binary2integer(binary)
        phenotype = _integer2phenotype(c, integer)
    end
    scalar = PhysicalFields.PhysicalScalar(phenotype, c.physical_units)
    return scalar
end # decode

# Methods to encode a parameter.

function _phenotype2integer(c::Chromosome, phenotype::Float64)::Int64
    fraction = (phenotype - c.parameter_min) / (c.parameter_max - c.parameter_min)
    integer  = Int(round(fraction * c.expressions))
    if integer < 0
        integer = 0
    end
    return convert(Int64, integer)
end # _phenotype2integer

function _integer2binary(c::Chromosome, integer::Int64)::Vector{Bool}
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

function encode!(c::Chromosome, scalar::PhysicalFields.PhysicalScalar)
    # verify input
    if c.parameter_units ≠ scalar.units
        msg = "The scalar's units are not equal to the chromosome's units."
        error(msg)
    end

    if c.genes == 0
        return nothing
    end

    phenotype = convert(Float64, PhysicalFields.get(scalar))
    if phenotype ≥ c.parameter_min && phenotype ≤ c.parameter_max
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
        msg = "The scalar must lie within [c.parameter_min, c.parameter_max]."
        error(msg)
    end
    return nothing
end # encode

# Methods for storing and retrieving a Chromosome to and from a file.

StructTypes.StructType(::Type{Chromosome}) = StructTypes.Struct()

"""
```
toFile(chromosome::GeneticAlgorithms.Chromosome, json_stream::IOStream)
```
writes a data structure `chromosome` to the IOStream `json_stream.`\n
For example, consider the code fragment:
```
json_stream = PhysicalFields.openJSONWriter(<my_dir_path>::String, <my_file_name>::String)\n
...\n
GeneticAlgorithms.toFile(chromosome::GeneticAlgorithms.Chromosome, json_stream::IOStream)\n
...\n
PhysicalFields.closeJSONStream(json_stream::IOStream)
```
where `<my_dir_path>` is the path to your working directory wherein the file\n
`<my_file_name>` that is to be written to either exists or will be created,\n
and which must have a .json extension.
"""
function toFile(chromosome::Chromosome, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, chromosome)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

"""
```
chromosome = fromFile(::GeneticAlgorithms.Chromosome, json_stream::IOStream)
```
reads an instance of type `Chromosome` from the IOStream `json_stream.`\n
For example, consider the code fragment:
```
json_stream = PhysicalFields.openJSONReader(<my_dir_path>::String, <my_file_name>::String)\n
...\n
chromosome = GeneticAlgorithms.fromFile(::GeneticAlgorithms.Chromosome, json_stream::IOStream)\n
...\n
PhysicalFields.closeJSONStream(json_stream::IOStream)\n
```
which returns a `chromosome,` an object of type `GeneticAlgorithms.Chromosome.`\n
Here `<my_dir_path>` is the path to your working directory wherein the file\n
to be read from, i.e., `<my_file_name>,` must exist, and which is to have a\n
.json extension.
"""
function fromFile(::Type{Chromosome}, json_stream::IOStream)::Chromosome
    if isopen(json_stream)
        chromosome = JSON3.read(readline(json_stream), Chromosome)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return chromosome
end
