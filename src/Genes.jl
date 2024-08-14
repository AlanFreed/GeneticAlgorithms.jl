# Admissible gene expressions include:

"""
```
const dominant
```

This constant is used to designate a dominant gene expression, e.g., in the constructor

```
my_dominant_gene = Gene(dominant)
```

"""
const dominant  = true    # Gene expression whenever a gene is dominant.

"""
```
const recessive
```

This constant is used to designate a recessive gene expression, e.g., in the constructor

```
my_recessive_gene = Gene(recessive)
```

"""
const recessive = false   # Gene expression whenever a gene is recessive.

"""
Genes are the lowest level containers of genetic information in a genetic algorithm.  Here haploid genes are considered, which have two expressions: dominant (assigned true) and recessive (assigned false).

The author considered diploid genes (which have three expressions: dominant, recessive, and dominant-recessive) early in his work with genetic algorithms, but found there to be no real advantage over using the simpler haploid gene.

A gene is a type whose datum is a gene expression.

```
struct Gene
    expression::PhysicalFields.MBoolean
end
```

where the `expression` is a mutable boolean, i.e., it can be changed/mutated.

Constructors

```
gene = Gene()
```
> Creates a `gene` with a random expression.

```
gene = Gene(dominant)
```
> Creates a `gene` with a `dominant` expression.

```
gene = Gene(recessive)
```
> Creates a `gene` with a `recessive` expression.

Operators

`==` and `≠`

Methods

```
e = get(g)
```
> Returns the expression `e` held by gene `g`.

```
set!(g, expression)
```
> Assigns a gene `expression` to the field `g.expression`.

```
c = copy(g)
```
> Returns a copy `c` of gene `g`.

```
s = toBinaryString(g)
```
> Returns a string representation `s` for gene `g` as either a "0" (recessive) or a "1" (dominant).

```
b = isDominant(g)
```
> Returns `b = true` if `g.expression == dominant`.

```
b = isRecessive(g)
```
> Returns `b = true` if `g.expression == recessive`.

```
mutate!(g, probability)
```
> Performs a random flip in the field `g.expression` at a specified `probability`, i.e., from `dominant` to `recessive`, or vice versa.

Persistence

```
toFile(g, json_stream)
```
> Writes a gene `g` to a JSON file `json_stream`.

```
g = fromFile(::Gene, json_stream)
```
> Reads a gene `g` from a JSON file `json_stream`.

Consider the following code fragments:

1) To open a file.

```
json_stream = PhysicalFields.openJSONWriter(<my_dir_path>, <my_file_name.json>)
```
> This opens a `json_stream` for `<my_file_name.json>` located in `<my_dir_path>`.

2) To write to a file.

```
toFile(gene, json_stream)
```
> This writes a `gene` to the `json_stream`.

3) To read from a file.

```
gene = GeneticAlgorithms.fromFile(::Gene, json_stream)
```
> This reads in a `gene` of type `Gene` from the `json_stream`.

4) And to close a file.

```
PhysicalFields.closeJSONStream(json_stream)
```
This flushes the buffer and closes the `json_stream`.
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

function fromFile(::Type{Gene}, json_stream::IOStream)::Gene
    if isopen(json_stream)
        gene = JSON3.read(readline(json_stream), Gene)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return gene
end
