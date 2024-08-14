# Admissible gene expressions include:

"""
```julia
const dominant
```
This constant is used to designate a dominant gene expression, e.g.,\
in the constructor\
```julia
my_dominant_gene = Gene(dominant)
```
"""
const dominant  = true    # Gene expression whenever a gene is dominant.

"""
```julia
const recessive
```
This constant is used to designate a recessive gene expression, e.g.,\
in the constructor\
```julia
my_recessive_gene = Gene(recessive)
```
"""
const recessive = false   # Gene expression whenever a gene is recessive.

"""
Genes are the lowest level containers of genetic information.  Here haploid\
genes are considered, which have two expressions: dominant (assigned true) and\
recessive (assigned false).\\
The author considered diploid genes (which have three expressions: dominant,\
recessive, and dominant-recessive) early in his work with genetic algorithms,\
but found there to be no advantage over using the simpler haploid gene.\\
A gene is a type whose datum is a gene expression.\
```julia
struct Gene
    expression::PhysicalFields.MBoolean
end
```
where the expression is a mutable boolean, i.e., it can be changed/mutated.\\

Constructors\
```julia
gene = Gene()  # Creates a `gene` with a random expression.
```
```julia
gene = Gene(dominant)  # Creates a `gene` with a dominant expression.
```
```julia
gene = Gene(recessive)  # Creates a `gene` with a recessive expression.
```
\\
Operators\

`==` and `≠`\\

Methods\

```julia
e = get(g)  # Returns expression `e` held by gene `g.`
```
```julia
set!(g, expression)  # Assigns gene `expression` to field `g.expression.`
```
```julia
c = copy(g)  # Returns a copy `c` of gene `g.`
```
```julia
s = toBinaryString(g)  # Returns a string representation `s` for gene `g.`
```
```julia
b = isDominant(g)  # Returns `b = true` if `g.expression == dominant.`
```
```julia
b = isRecessive(g)  # Returns `b = true` if `g.expression == recessive.`
```
```julia
mutate!(g, probability)  # A random flip in `g.expression` at a `probability.`
```
```julia
toFile(g, json_stream)  # Writes gene `g` to a JSON file.
```
```julia
g = fromFile(::Gene, json_stream)  # Reads gene `g` from a JSON file.
"""
struct Gene
    expression::PhysicalFields.MBoolean

    # constructors

    function Gene(expression::PhysicalFields.MBoolean)
        new(expression)::Gene
    end

    function Gene(expression::Bool)
        gene_expression = PhysicalFields.MBoolean(expression)
        new(gene_expression)::Gene
    end

    function Gene()
        value = rand()    # Returns a random real from the interval [0, 1).
        if value < 0.5
            expression = PhysicalFields.MBoolean(recessive)
        else
            expression = PhysicalFields.MBoolean(dominant)
        end
        new(expression)::Gene
    end
end # Gene

# operators

function Base.:(==)(gL::Gene, gR::Gene)::Bool
    if gL.expression == gR.expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(g::Gene, expression::MBoolean)::Bool
    if g.expression == expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(expression::MBoolean, g::Gene)::Bool
    if expression == g.expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(g::Gene, expression::Bool)::Bool
    if g.expression == expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(expression::Bool, g::Gene)::Bool
    if expression == g.expression
        return true
    else
        return false
    end
end # ==

function Base.:≠(gL::Gene, gR::Gene)::Bool
    if gL.expression ≠ gR.expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(g::Gene, expression::MBoolean)::Bool
    if g.expression ≠ expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(expression::MBoolean, g::Gene)::Bool
    if expression ≠ g.expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(g::Gene, expression::Bool)::Bool
    if g.expression ≠ expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(expression::Bool, g::Gene)::Bool
    if expression ≠ g.expression
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(get)(g::Gene)::Bool
    return get(g.expression)
end # get

function set!(g::Gene, expression::Bool)
    set!(g.expression, expression)
    return nothing
end # set!

function Base.:(copy)(g::Gene)::Gene
    expression = copy(g.expression)
    return Gene(expression)
end # copy

function toBinaryString(g::Gene)::String
    if g == recessive
        s = "0"
    else
        s = "1"
    end
    return s
end # toBinaryString

function isDominant(g::Gene)::Bool
    if g == dominant
        return true
    else
        return false
    end
end # isDominant

function isRecessive(g::Gene)::Bool
    if g == recessive
        return true
    else
        return false
    end
end # isRecessive

"""
mutate! addresses the occurrence of a gene expression changing from dominant\
to recessive, or vice versa, with a chance of change at `probability_mutation.`
"""
function mutate!(g::Gene, probability_mutation::Float64)
    if probability_mutation < 0.0 || probability_mutation ≥ 1.0
        msg = "A probability of mutation must belong to unit interval [0, 1)."
        error(msg)
    end
    if probability_mutation > rand()
        if isDominant(g)
            set!(g, recessive)
        else
            set!(g, dominant)
        end
    end
    return nothing
end # mutate!

function mutate!(g::Gene, probability_mutation::Real)
    probability = convert(Float64, probability_mutation)
    mutate!(g, probability)
    return nothing
end

# Methods for storing and retrieving a Gene to and from a file.

StructTypes.StructType(::Type{Gene}) = StructTypes.Struct()

"""
```julia
toFile(gene::GeneticAlgorithms.Gene, json_stream::IOStream)
```
writes a data structure `gene` to the IOStream `json_stream.`\\
For example, consider the code fragment:
```julia
json_stream = PhysicalFields.openJSONWriter(<my_dir_path>::String, <my_file_name>::String)\
...\
GeneticAlgorithms.toFile(gene::GeneticAlgorithms.Gene, json_stream::IOStream)\
...\
PhysicalFields.closeJSONStream(json_stream::IOStream)
```
where `<my_dir_path>` is the path to your working directory wherein the file\
`<my_file_name>` that is to be written to either exists or will be created,\
and which must have a .json extension.
"""
function toFile(gene::Gene, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, gene)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

"""
```julia
gene = fromFile(::GeneticAlgorithms.Gene, json_stream::IOStream)
```
reads an instance of type `Gene` from the IOStream `json_stream.`\
For example, consider the code fragment:
```julia
json_stream = PhysicalFields.openJSONReader(<my_dir_path>::String, <my_file_name>::String)\
...\
gene = GeneticAlgorithms.fromFile(::GeneticAlgorithms.Gene, json_stream::IOStream)\
...\
PhysicalFields.closeJSONStream(json_stream::IOStream)
```
which returns a `gene,` an object of type `GeneticAlgorithms.Gene.`\
Here `<my_dir_path>` is the path to your working directory wherein the file\
to be read from, i.e., `<my_file_name>,` must exist, and which is to have a\
.json extension.
"""
function fromFile(::Type{Gene}, json_stream::IOStream)::Gene
    if isopen(json_stream)
        gene = JSON3.read(readline(json_stream), Gene)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return gene
end
