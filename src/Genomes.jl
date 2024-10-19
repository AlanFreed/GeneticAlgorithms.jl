"""
A genome is a genetic container of chromosomes.  A genome represents all parameters in a model being parameterized in terms of their parameters' genetic counterparts.

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# Genome

```julia
struct Genome
    parameters_constrained::Vector{Tuple{Int,Int}}
    genes::Int
    chromosomes::Int
    genotypes::Vector{Chromosome}
end
```
where
1. *parameters_constrained* provides parameter constraints as tuple pairs, e.g., tuple (l, r) imposes a constraint of θ[l] < θ[r].
2. *genes* specifies the number of genes that comprise a genome.
3. *chromosomes* specifies the number of chromosomes that comprise a genome.
4. *genotypes* provides the genetic material, viz., the chromosomes comprising a genome.

## Constructor

```julia
g = Genome(parameters_min::Vector{PF.PhysicalScalar},
           parameters_max::Vector{PF.PhysicalScalar},
           parameters_constrained::Vector{Tuple{Int,Int}}, 
           significant_figures::Int)
```
where
1. *parameters_min* is a vector containing the minimum values for each phenotype (parameter in the model) being parameterized, i.e., θₘᵢₙ[i].
2. *parameters_max* is a vector containing the maximum values for each phenotype in the model being parameterized, i.e., θₘₐₓ[i].
3. *parameters_constrained* is a vector containing constraints between parameters, if any exit. If not, the user should send `parameters_constrained = Vector{Tuple}(undef, 0)` to this constructor.
4. *significant_figures* specifies the number of digits in accuracy sought by a solution for the model's parameters.

> **NOTE**: if `parameters_min[i] ≈ parameters_max[i]` then the phenotype being represented by the chromosome at index i will be taken to be constant, fixed to the value of this collapsed range.

## Operators

`==` and `≠`

## Methods

```julia
chromosome = getindex(genome::Genome, index::Int)
```
> returns a chromosome from a genome at location index, or `chromosome = genome[index]`.

```julia
setindex!(genome::Genome, chromosome::Chromosome, index::Int)
```
> assigns a chromosome to a genome at location index, or `genome[index] = chromosome`.

```julia
cc = copy(genome::Genome)
```
> returns a copy cc of the genome.

```julia
str = toBinaryString(genome::Genome)
```
> returns a string str describing the genome written in a binary format.

```julia
str = toString(genome::Genome)
```
> returns a string str representing the phenotypes expressed by a genome. 

```julia
mutate!(genome::Genome, probability_mutation::Float64)
```
> performs a random flip in gene expression (i.e., dominant to recessive or vice versa) at a specified probability of mutation applied to each gene in the genome.

```julia
C = crossover(A::Genome, B::Genome, probability_mutation::Float64, probability_crossover::Float64)
```
> Performs a crossover or conception between two genomes coming from *parents* A and B with a possibility for individual gene mutations and with a likelyhood for crossover (chromosome splitting or sharing between parents).  The result is a *child* genome C.

```julia
θ = decode(genome::Genome)
```
> returns a vector of phenotypes (the parameters θ of a model: a vector of physical scalars, possibly with different physical units, i.e., an instance of type `Vector{PF.PhysicalScalar}`) held by the genome.

```julia
encode!(genome::Genome, θ::Vector{PF.PhysicalScalar})
```
> assigns phenotypes (the parameters θ of a model) to a genome.

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

To write or read an instance of type *Genome* to or from a JSON file, call
```julia
toFile(genome, json_stream)
```
> which writes a genome of type *Genome* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
genome = fromFile(Genome, json_stream)
```
> reads a genome of type *Genome* from the JSON file attached to json_stream.
"""
struct Genome
    # Field that helps describe the parameters.
    parameters_constrained::Vector{Tuple{Int,Int}}
    # Fields that describe the genome.
    genes::Int
    chromosomes::Int
    genotypes::Vector{Chromosome}

    # constructors

    function Genome(parameters_min::Vector{PF.PhysicalScalar},
                    parameters_max::Vector{PF.PhysicalScalar}, 
                    parameters_constrained::Vector{Tuple{Int,Int}}, 
                    significant_figures::Int)
        # verify input
        if length(parameters_min) == length(parameters_max)
            chromosomes = length(parameters_min)
        else
            msg = "Vectors parameters_min and parameters_max "
            msg = string(msg, "must have the same length.")
            throw(DimensionMismatch, msg)
        end
        for i in 1:length(parameters_min)
            p_min = parameters_min[i]
            p_max = parameters_max[i]
            if p_min.units ≠ p_max.units
                msg = string("parameters_min[", i, "] and ")
                msg = string(msg, "parameters_max[", i, "]\n")
                msg = string(msg, "must have the same physical units.")
                error(msg)
            end
        end

        # create the genome with randomly assigned genes
        genotypes = Vector{Chromosome}(undef, chromosomes)
        for i in 1:chromosomes
            genotypes[i] = Chromosome(parameters_min[i],
                                      parameters_max[i], 
                                      significant_figures)
        end

        genes = 0
        for i in 1:chromosomes
            chromosome = genotypes[i]
            genes = genes + chromosome.genes
        end

        new(parameters_constrained, genes, chromosomes, genotypes)::Genome
    end
    
    function Genome(parameters_constrained::Vector{Tuple{Int,Int}},
                    genes::Int,
                    chromosomes::Int, 
                    genotypes::Vector{Chromosome})

        new(parameters_constrained, genes, chromosomes, genotypes)::Genome
    end
end # Genome

# operators

function Base.:(==)(genome_left::Genome, genome_right::Genome)::Bool
    if genome_left.genes ≠ genome_right.genes
        return false
    end
    if genome_left.chromosomes ≠ genome_right.chromosomes
        return false
    end
    for i in 1:length(genome_left.parameters_constrained)
        if (genome_left.parameters_constrained[i] ≠ 
            genome_right.parameters_constrained[i])
            return false
        end
    end
    for i in 1:genome_left.chromosomes
        if genome_left.genotypes[i] ≠ genome_right.genotypes[i]
            return false
        end
    end
    return true
end # ==

function Base.:≠(genome_left::Genome, genome_right::Genome)::Bool
    return !(genome_left == genome_right)
end # ≠

# methods

function Base.:(getindex)(genome::Genome, index::Int)::Chromosome
    if index ≥ 1 && index ≤ genome.chromosomes
        chromosome = genome.genotypes[index]
    else
        msg = "Admissible chromosome indices are ∈ [1…"
        msg = string(msg, genome.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return copy(chromosome)
end # getindex

function Base.:(setindex!)(genome::Genome, chromosome::Chromosome, index::Int)
    if index ≥ 1 && index ≤ genome.chromosomes
        old_chromosome = genome.genotypes[index]
        if old_chromosome.parameter_max.units ≠ chromosome.parameter_max.units
            error("Chromosome to be assigned has the wrong physical units.")
        end
        genome.genotypes[index] = copy(chromosome)
    else
        msg = "Admissible chromosome indices are ∈ [1…"
        msg = string(msg, genome.chromosomes, "].")
        throw(DimensionMismatch(msg))
    end
    return nothing
end # setindex!

function Base.:(copy)(genome::Genome)::Genome
    constraints = length(genome.parameters_constrained)
    parameters_constrained = Vector{Tuple{Int,Int}}(undef, constraints)
    for i in 1:constraints
        (left, right) = genome.parameters_constrained[i]
        parameters_constrained[i] = (copy(left), copy(right))
    end
    genes = copy(genome.genes)
    chromosomes = copy(genome.chromosomes)
    genotypes = Vector{Chromosome}(undef, chromosomes)
    for i in 1:chromosomes
        genotypes[i] = copy(genome.genotypes[i])
    end
    return Genome(parameters_constrained, genes, chromosomes, genotypes)
end # copy

function toBinaryString(genome::Genome)::String
    s = ""
    for i in 1:genome.chromosomes
        chromosome = genome.genotypes[i]
        if genome.chromosomes < 10
            s = string(s, i, ": ", toBinaryString(chromosome))
        elseif genome.chromosomes < 100
            if i < 10
                s = string(s, i, ":  ", toBinaryString(chromosome))
            else
                s = string(s, i, ": ", toBinaryString(chromosome))
            end
        else
            if i < 10
                s = string(s, i, ":   ", toBinaryString(chromosome))
            elseif i < 100
                s = string(s, i, ":  ", toBinaryString(chromosome))
            else
                s = string(s, i, ": ", toBinaryString(chromosome))
            end
        end
        if i < genome.chromosomes
            s = string(s, "\n")
        end
    end
    return s
end # toBinaryString

function toString(genome::Genome)::String
    s = ""
    for i in 1:genome.chromosomes
        chromosome = genome.genotypes[i]
        phenotype  = decode(chromosome)
        if genome.chromosomes < 10
            if phenotype.value > -0.0
                s = string(s, i, ":  ", PF.toString(phenotype))
            else
                s = string(s, i, ": ", PF.toString(phenotype))
            end
        elseif genome.chromosomes < 100
            if phenotype.value > -0.0
                if i < 10
                    s = string(s, i, ":    ", PF.toString(phenotype))
                elseif i < 100
                    s = string(s, i, ":   ", PF.toString(phenotype))
                else
                    s = string(s, i, ":  ", PF.toString(phenotype))
                end
            else
                if i < 10
                    s = string(s, i, ":   ", PF.toString(phenotype))
                elseif i < 100
                    s = string(s, i, ":  ", PF.toString(phenotype))
                else
                    s = string(s, i, ": ", PF.toString(phenotype))
                end
            end
        end
        if i < genome.chromosomes
            s = string(s, "\n")
        end
    end
    return s
end # toString

function mutate!(genome::Genome, probability_mutation::Float64)
    for i in 1:genome.chromosomes
        chromosome = genome.genotypes[i]
        mutate!(chromosome, probability_mutation)
        genome.genotypes[i] = chromosome
    end
    return nothing
end # mutate!

function crossover(parentA::Genome, parentB::Genome,
                   probability_mutation::Float64,
                   probability_crossover::Float64)::Genome

    if parentA.genes == parentB.genes
        if parentA.chromosomes == parentB.chromosomes
            child = copy(parentA)
            for i in 1:child.chromosomes
                child.genotypes[i] = crossover(parentA.genotypes[i],
                                               parentB.genotypes[i],
                                               probability_mutation,
                                               probability_crossover)
            end
        else
            msg = "Parents must have the same number of chromosomes."
            throw(DimensionMismatch(msg))
        end
    else
        msg = "Parents must have the same number of genes."
        throw(DimensionMismatch(msg))
    end
    return child
end # crossover

# Cannot return ArrayOfPhysicalScalars because units may vary between entries.
function decode(genome::Genome)::Vector{PF.PhysicalScalar}
    phenotypes = Vector{PF.PhysicalScalar}(undef, genome.chromosomes)
    for chromosome in 1:genome.chromosomes
        phenotypes[chromosome] = decode(genome.genotypes[chromosome])
    end
    return phenotypes
end # decode

function encode!(genome::Genome, phenotypes::Vector{PF.PhysicalScalar})
    # verify input
    if genome.chromosomes ≠ length(phenotypes)
        msg = "Number of chromosomes must equal number of phenotypes."
        throw(DimensionMismatch, msg)
    end

    for i in 1:genome.chromosomes
        chromosome = genome.genotypes[i]
        phenotype  = phenotypes[i]
        encode!(chromosome, phenotype)
        genome.genotypes[i] = chromosome
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

