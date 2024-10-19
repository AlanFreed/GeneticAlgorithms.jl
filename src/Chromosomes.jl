"""
A chromosome is a genetic container of genes. In this implementation of a genetic algorithm, each chromosome represents a parameter (an unknown scalar in some model to be parameterized). Genetic processes (mutation and crossover) adjust these chromosomes through evolutionary processes, and therefore, their associated parametric values change over the generations. Chromosomes are gray binary representations for a model's parameters θ.

Chromosomes are where genetics and optimization meet.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# Chromosome

```julia
struct Chromosome
    parameter_min::PF.PhysicalScalar
    parameter_max::PF.PhysicalScalar
    genes::Int
    expressions::Int64
    genotype::Vector{Gene}
end
```
where
1. *parameter_min* specifies the minimum scalar value that the phenotype (value of the parameter) can represent, i.e., θₘᵢₙ.
2. *parameter_max* specifies the maximum scalar value that the phenotype can represent, i.e., θₘₐₓ.
3. *genes* specifies the number of genes that comprise a chromosome.
4. *expressions* specifies the number of gene expressions a chromosome can represent.
5. *genotype* provides the genetic material, viz., the genes comprising a chromosome.

A parameter will be considered fixed if its `parameter_min ≈ parameter_max`.

## Constructor

The constructor most likely to be called by the user is
```julia
chromosome = Chromosome(parameter_min::PF.PhysicalScalar,
                        parameter_max::PF.PhysicalScalar,
                        significant_figures::Int)
```
> where `parameter_min` and `parameter_max` are described above, while `significant_figures` specifies the number of digits in accuracy sought by a solution for the model's parameters.

The general constructor called, e.g., by JSON3 is
```julia
chromosome = Chromosome(parameter_min::PF.PhysicalScalar, 
                        parameter_max::PF.PhysicalScalar,
                        genes::Int, 
                        expressions::Int64, 
                        genotype::Vector{Gene})
```

> **NOTE**: if `parameter_min ≈ parameter_max` then the phenotype being represented by this chromosome will be taken to be constant, fixed to the value of this collapsed range.

## Operators

`==` and `≠`

## Methods

```julia
gene = getindex(chromosome::Chromosome, index::Int)
```
> returns a gene from a chromosome located at an index, or `gene = chromosome[index]`.

```julia
setindex!(chromosome::Chromosome, gene::Gene, index::Int)
```
> assigns a gene to a chromosome  located at an index, or `chromosome[index] = gene`.

```julia
cc = copy(c::Chromosome)
```
> returns a copy cc of a chromosome.

```julia
str = toBinaryString(chromosome::Chromosome)
```
> returns a string str that describes a chromosome written in a binary format.

```julia
str = toString(chromosome::Chromosome)
```
> returns a string str for the phenotype represented by a chromosome. 

```julia
mutate!(chromosome::Chromosome, probability_mutation::Float64)
```
> performs a random flip in gene expression (i.e., dominant to recessive or vice versa) at a specified probability for mutation applied to each gene in a chromosome.

```julia
C = crossover(A::Chromosome, B::Chromosome, probability_mutation::Float64, probability_crossover::Float64)
```
> performs a crossover at conception between two chromosomes coming from parents A and B with possibilities for individual gene mutation and crossover (chromosome splitting). The result is a child chromosome C.

```julia
θ = decode(chromosome::Chromosome)
```
> returns a phenotype (the parameter θ, which is a PF.PhysicalScalar) held in the gene expression of a chromosome.

```julia
encode!(chromosome::Chromosome, θ::PF.PhysicalScalar)
```
> assigns a phenotype θ (value of a parameter) to a chromosome.

### Persistence

To open or close an IOStream attached to a JSON file, call
```julia
json_stream = PF.openJSONWriter(<my_dir_path>, <my_file_name.json>)
```
> which opens a `json_strea`m of type *IOStream* for a file `<my_file_name.json>` located in directory `<my_dir_path>`, both of which are strings, while
```julia
PF.closeJSONStream(json_stream)
```
> flushes the buffer and closes this json_stream.

To write or read an instance of type *Chromosome* to or from a JSON file, call
```julia
toFile(chromosome, json_stream)
```
> which writes a chromosome of type *Chromosome* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
chromosome = fromFile(Chromosome, json_stream)
```
> reads a chromosome of type *Chromosome* from the JSON file attached to `json_stream`.
"""
struct Chromosome
    # Fields that describe a parameter.
    parameter_min::PF.PhysicalScalar
    parameter_max::PF.PhysicalScalar
    # Fields that describe a chromosome.
    genes::Int
    expressions::Int64
    genotype::Vector{Gene}

    # constructor

    function Chromosome(parameter_min::PF.PhysicalScalar, 
                        parameter_max::PF.PhysicalScalar,
                        significant_figures::Int)
        if parameter_min.units ≠ parameter_max.units
            error("A parameter's bounds must have the same physical units.")
        end
        
        if parameter_min ≈ parameter_max
            # Handle the case of a fixed parameter.
            parameter_min = parameter_max
            genes         = 0
            expressions   = Int64(1)
            genotype      = Vector{Gene}(undef, 0)
        else
            # Verify the input.
            if parameter_min > parameter_max
                msg = "Cannot create a chromosome unless parameter_max "
                msg = string(msg, "≥ parameter_min.")
                error(msg)
            end

            # Number of decades that span [parameter_min, parameter_max].
            if parameter_min.value > 0.0
                log_decades = log10(parameter_max.value/parameter_min.value)
            elseif parameter_max.value < 0.0
                log_decades = log10(parameter_min.value/parameter_max.value)
            elseif parameter_min.value ≈ 0.0
                log_decades = log10(parameter_max.value)
            elseif parameter_max.value ≈ 0.0
                log_decades = log10(-parameter_min.value)
            else # parameter_min < 0.0 and parameter_max > 0.0
                log_decades = log10(-parameter_min.value*parameter_max.value)
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
            genes = genes_per_decade * decades
            if genes < 63
                expressions = Int64(2^genes - 1)
            else
                genes = 63
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

        new(parameter_min, parameter_max,
            genes, expressions, genotype)::Chromosome
    end

    function Chromosome(parameter_min::PF.PhysicalScalar, 
                        parameter_max::PF.PhysicalScalar,
                        genes::Int, 
                        expressions::Int64, 
                        genotype::Vector{Gene})

        new(parameter_min, parameter_max,
            genes, expressions, genotype)::Chromosome
    end
end # Chromosome

# operators

function Base.:(==)(chromosome_left::Chromosome,
                    chromosome_right::Chromosome)::Bool
    if (chromosome_left.parameter_max.units ≠ 
        chromosome_right.parameter_max.units || 
        chromosome_left.genes ≠ chromosome_right.genes)
        return false
    end
    if (chromosome_left.genes == 0 && 
        chromosome_left.parameter_min ≈ chromosome_right.parameter_min)
        return true
    end
    for i in 1:chromosome_left.genes
        if chromosome_left.genotype[i] ≠ chromosome_right.genotype[i]
            return false
        end
    end
    return true
end # ==

function Base.:≠(chromosome_left::Chromosome,
                 chromosome_right::Chromosome)::Bool
    return !(chromosome_left == chromosome_right)
end # ≠

# methods

function Base.:(getindex)(chromosome::Chromosome, index::Int)::Gene
    if chromosome.genes == 0
        return nothing
    elseif index ≥ 1 && index ≤ chromosome.genes
        gene = chromosome.genotype[index]
    else
        msg = "Admissible gene indices are ∈ [1…"
        msg = string(msg, chromosome.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return copy(gene)
end # getindex

function Base.:(setindex!)(chromosome::Chromosome, gene::Gene, index::Int)
    if chromosome.genes == 0
        return nothing
    elseif index ≥ 1 && index ≤ chromosome.genes
        set!(chromosome.genotype[index], get(gene))
    else
        msg = "Admissible gene indices are ∈ [1…"
        msg = string(msg, chromosome.genes, "].")
        throw(DimensionMismatch(msg))
    end
    return nothing
end # setindex!

function Base.:(copy)(chromosome::Chromosome)::Chromosome
    parameter_min   = PF.copy(chromosome.parameter_min)
    parameter_max   = PF.copy(chromosome.parameter_max)
    genes           = copy(chromosome.genes)
    expressions     = copy(chromosome.expressions)
    genotype        = Vector{Gene}(undef, genes)
    for i in 1:genes
        genotype[i] = copy(chromosome.genotype[i])
    end
    return Chromosome(parameter_min, parameter_max, 
                      genes, expressions, genotype)
end # copy

function toBinaryString(chromosome::Chromosome)::String
    s = ""
    for i in 1:chromosome.genes
        gene = chromosome.genotype[i]
        s = string(s, toBinaryString(gene))
    end
    return s
end # toBinaryString

function toString(chromosome::Chromosome)::String
    return PF.toString(decode(chromosome))
end # toString

function mutate!(chromosome::Chromosome, probability_mutation::Float64)
    for i in 1:chromosome.genes
        gene = chromosome.genotype[i]
        mutate!(gene, probability_mutation)
        chromosome.genotype[i] = gene
    end
    return nothing
end # mutate!

function crossover(parentA::Chromosome, parentB::Chromosome,
                   probability_mutation::Float64,
                   probability_crossover::Float64)::Chromosome

    if probability_crossover < 0.0 || probability_crossover ≥ 1.0
        error("Probability of crossover must belong to unit interval [0, 1).")
    end

    if parentA.genes == parentB.genes
        if rand() < 0.5
            child = copy(parentA)
            # Left side of chromosome splice belongs to parent A.
            if probability_crossover > rand()
                if child.genes > 3
                    xover = rand(2:child.genes-1)
                    for i in xover:child.genes
                        child.genotype[i] = parentB.genotype[i]
                    end
                elseif child.genes == 3
                    if rand() < 0.1
                        child.genotype[2] = parentB.genotype[2]
                    end
                    child.genotype[3] = parentB.genotype[3]
                elseif child.genes == 2
                    child.genotype[2] = parentB.genotype[2]
                else
                    # do nothing
                end
            end
        else
            child = copy(parentB)
            # Left side of chromosome splice belongs to parent B.
            if probability_crossover > rand()
                if child.genes > 3
                    xover = rand(2:child.genes-1)
                    for i in xover:child.genes
                        child.genotype[i] = parentA.genotype[i]
                    end
                elseif child.genes == 3
                    if rand() < 0.5
                        child.genotype[2] = parentA.genotype[2]
                    end
                    child.genotype[3] = parentA.genotype[3]
                elseif child.genes == 2
                    child.genotype[2] = parentA.genotype[2]
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
This algorithm is for mapping between binary and gray codes.  It considers the
most significant bit (MSB) is at the [Low] position of the code, while the 
least significant bit (LSB) associates with the [High] position of the code, 
in other words, e.g.,
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
    integer = Int64(0)
    power   = Int64(1)
    for i in bits:-1:1
        if binary[i] == dominant
            integer = integer + power
        end
        power = 2power
    end
    return integer
end # _binary2integer

function _integer2phenotype(chromosome::Chromosome,
                            integer::Int64)::PF.PhysicalScalar
    phenotype = (chromosome.parameter_min + 
                 Float64(integer / chromosome.expressions) *
                 (chromosome.parameter_max - chromosome.parameter_min))
    if phenotype < chromosome.parameter_min
        phenotype = chromosome.parameter_min
    end
    if phenotype > chromosome.parameter_max
        phenotype = chromosome.parameter_max
    end
    return phenotype
end # _integer2phenotype

function decode(chromosome::Chromosome)::PF.PhysicalScalar
    if chromosome.genes == 0
        phenotype = copy(chromosome.parameter_min)
    else
        gray = Vector{Bool}(undef, chromosome.genes)
        for i in 1:chromosome.genes
            if isDominant(chromosome.genotype[i])
                gray[i] = dominant
            else
                gray[i] = recessive
            end
        end
        binary    = _gray2binary(gray)
        integer   = _binary2integer(binary)
        phenotype = _integer2phenotype(chromosome, integer)
    end
    return phenotype
end # decode

# Methods to encode a parameter.

function _phenotype2integer(chromosome::Chromosome,
                            phenotype::PF.PhysicalScalar)::Int64
    fraction = PF.get((phenotype - chromosome.parameter_min) /
                      (chromosome.parameter_max - chromosome.parameter_min))
    integer  = Int64(round(fraction*chromosome.expressions))
    if integer < 0
        integer = Int64(0)
    end
    return integer
end # _phenotype2integer

function _integer2binary(chromosome::Chromosome, integer::Int64)::Vector{Bool}
    binary    = Vector{Bool}(undef, chromosome.genes)
    atinteger = integer
    gene      = chromosome.genes
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

function encode!(chromosome::Chromosome, phenotype::PF.PhysicalScalar)
    # verify input
    if chromosome.parameter_max.units ≠ phenotype.units
        error("Phenotype's units are not equal to chromosome's units.")
    end

    # handle the constant parameter
    if chromosome.genes == 0
        return nothing
    end
    
    # encode
    if (phenotype ≥ chromosome.parameter_min && 
        phenotype ≤ chromosome.parameter_max)
        integer = _phenotype2integer(chromosome, phenotype)
        binary  = _integer2binary(chromosome, integer)
        gray    = _binary2gray(binary)
        for i in 1:chromosome.genes
            if gray[i] == dominant
                chromosome.genotype[i] = Gene(dominant)
            else
                chromosome.genotype[i] = Gene(recessive)
            end
        end
    else
        msg = "The phenotype's value must lie within ["
        msg = string(msg, chromosome.parameter_min.value, ", ")
        msg = string(msg, chromosome.parameter_max.value, "] ")
        msg = string(msg, PF.toString(phenotype.units), ".")
        error(msg)
    end
    return nothing
end # encode

# Methods for storing and retrieving a Chromosome to and from a file.

StructTypes.StructType(::Type{Chromosome}) = StructTypes.Struct()

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

function fromFile(::Type{Chromosome}, json_stream::IOStream)::Chromosome
    if isopen(json_stream)
        chromosome = JSON3.read(readline(json_stream), Chromosome)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return chromosome
end

