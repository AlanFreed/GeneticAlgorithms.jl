# Admissible gene expressions include:

"""
```julia
const dominant
```
This constant is used to designate a dominant gene expression, e.g., in the constructor
```julia
my_dominant_gene = Gene(dominant)
```
"""
const dominant  = true

"""
```julia
const recessive
```
This constant is used to designate a recessive gene expression, e.g., in the constructor
```julia
my_recessive_gene = Gene(recessive)
```
"""
const recessive = false

# Type Gene

"""
Genes are the lowest level containers of genetic information in a genetic algorithm. Here haploid genes are considered. They have two expressions: dominant (assigned as `true`) and recessive (assigned as `false`).

To simplify this help, we use the alias
```julia
import
    PhysicalFields as PF
```

# Gene

```julia
struct Gene
    expression::PF.MBoolean
end
```
where the *expression* is a mutable boolean, i.e., it can be changed/mutated.

## Constructors

```julia
gene = Gene()
```
> creates a gene with a random expression.

```julia
gene = Gene(dominant)
```
> creates a gene with a dominant expression.

```julia
gene = Gene(recessive)
```
> creates a gene with a recessive expression.

## Operators

`==` and `≠`

## Methods

```julia
expression = get(gene::Gene)
```
> returns the expression held by a gene, i.e., gene.expression.

```julia
set!(gene::Gene, expression::Bool)
```
> assigns a gene expression to field gene.expression.

```julia
cc = copy(gene::Gene)
```
> returns a copy cc of a gene.

```julia
str = toBinaryString(gene::Gene)
```
> returns a string representation str for a gene as either a "0" (recessive) or a "1" (dominant).

```julia
boolean = isDominant(gene::Gene)
```
> returns `boolean = true` if `gene.expression == dominant`.

```julia
boolean = isRecessive(gene::Gene)
```
> returns `boolean = true` if `gene.expression == recessive`.

```julia
mutate!(gene::Gene, probability::Float64)
```
> Performs a random flip in the field gene.expression at a specified  probability, i.e., from dominant to recessive, or vice versa.

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

To write or read an instance of type *Gene* to or from a JSON file, call
```julia
toFile(gene, json_stream)
```
> which writes a gene of type *Gene* to the JSON file attached to a `json_stream` of type *IOStream*, while
```julia
gene = fromFile(Gene, json_stream)
```
> reads a gene of type *Gene* from the JSON file attached to json_stream.
"""
struct Gene
    expression::PF.MBoolean

    # constructors

    function Gene(expression::PF.MBoolean)
        new(expression)::Gene
    end

    function Gene(expression::Bool)
        gene_expression = PF.MBoolean(expression)
        new(gene_expression)::Gene
    end

    function Gene()
        value = rand()  # Returns a random real from the interval [0, 1).
        if value < 0.5
            expression = PF.MBoolean(recessive)
        else
            expression = PF.MBoolean(dominant)
        end
        new(expression)::Gene
    end
end # Gene

# operators

function Base.:(==)(gene_left::Gene, gene_right::Gene)::Bool
    if gene_left.expression == gene_right.expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(gene::Gene, expression::MBoolean)::Bool
    if gene.expression == expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(expression::MBoolean, gene::Gene)::Bool
    if expression == gene.expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(gene::Gene, expression::Bool)::Bool
    if gene.expression == expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(expression::Bool, gene::Gene)::Bool
    if expression == gene.expression
        return true
    else
        return false
    end
end # ==

function Base.:≠(gene_left::Gene, gene_right::Gene)::Bool
    return !(gene_left == gene_right)
end # ≠

function Base.:≠(gene::Gene, expression::MBoolean)::Bool
    return !(gene == expression)
end # ≠

function Base.:≠(expression::MBoolean, gene::Gene)::Bool
    return !(expression == gene)
end # ≠

function Base.:≠(gene::Gene, expression::Bool)::Bool
    return !(gene == expression)
end # ≠

function Base.:≠(expression::Bool, gene::Gene)::Bool
    return !(expression == gene)
end # ≠

# methods

function Base.:(get)(gene::Gene)::Bool
    return PF.get(gene.expression)
end # get

function set!(gene::Gene, expression::Bool)
    PF.set!(gene.expression, expression)
    return nothing
end # set!

function Base.:(copy)(gene::Gene)::Gene
    expression = copy(gene.expression)
    return Gene(expression)
end # copy

function toBinaryString(gene::Gene)::String
    if gene == recessive
        s = "0"
    else
        s = "1"
    end
    return s
end # toBinaryString

function isDominant(gene::Gene)::Bool
    if gene == dominant
        return true
    else
        return false
    end
end # isDominant

function isRecessive(gene::Gene)::Bool
    if gene == recessive
        return true
    else
        return false
    end
end # isRecessive

function mutate!(gene::Gene, probability_mutation::Float64)
    if probability_mutation < 0.0 || probability_mutation ≥ 1.0
        msg = "Probability for mutation must belong to unit interval [0, 1)."
        error(msg)
    end
    if probability_mutation > rand()
        if isDominant(gene)
            set!(gene, recessive)
        else
            set!(gene, dominant)
        end
    end
    return nothing
end # mutate!

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

