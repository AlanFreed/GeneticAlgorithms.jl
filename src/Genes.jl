# Admissible gene expressions include:

const dominant  = true    # Gene expression whenever a gene is dominant.
const recessive = false   # Gene expression whenever a gene is recessive.

"""
Genes are the lowest level containers of genetic information.  Here they are
considered to contain haploid gene expressions, i.e., dominant (assigned true)
and recessive (assigned false).

A gene is a type whose datum is a gene expression.

struct Gene
    expression::PhysicalFields.MBoolean
end

Constructors

    g = Gene()                  creates a gene 'g' with a random expression
    g = Gene(dominant)          creates a gene 'g' with dominant expression
    g = Gene(recessive)         creates a gene 'g' with recessive expression

Operators

    ==, ≠

Methods

    e = get(g)                  returns expression 'e' held by gene 'g'
    set!(g, expression)         assigns a gene `expression` to g.expression
    c = copy(g)                 returns a copy 'c' of gene 'g'
    s = toBinaryString(g)       returns a string representation 's' of gene 'g'
    b = isDominant(g)           returns true if g.expression == dominant
    b = isRecessive(g)          returns true if g.expression == recessive
    mutate!(g, probability)     random flip in gene expression at 'probability'
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
mutate! addresses the occurrence of a gene's expression changing from dominant
to recessive, or vice versa, with a chance of change at probability_mutation.
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
Method:\n
    toFile(g::GeneticAlgorithms.Gene, json_stream::IOStream)\n
Writes data structure `g` to the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONWriter(<my_dir_path>, <my_file_name>)\n
    ...\n
    GeneticAlgorithm.toFile(g::Gene, json_stream::IOStream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream::IOStream)\n
where <my_dir_path> is the path to your working directory wherein the file
<my_file_name> that is to be written to either exists or will be created, and
which must have a .json extension.
"""
function toFile(g::Gene, json_stream::IOStream)
    if isopen(json_stream)
        JSON3.write(json_stream, g)
        write(json_stream, '\n')
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    flush(json_stream)
    return nothing
end

"""
Method:\n
    fromFile(g::GeneticAlgorithms.Gene, json_stream::IOStream)\n
Reads a Gene from the IOStream `json_stream.`\n
For example, consider the code fragment:\n
    json_stream = PhysicalFields.openJSONReader(<my_dir_path>, <my_file_name>)\n
    ...\n
    g = GeneticAlgorithms.fromFile(GeneticAlgorithms.Gene, json_stream)\n
    ...\n
    PhysicalFields.closeJSONStream(json_stream)\n
that returns `g,` which is an object of type GeneticAlgorithms.Gene. Here
<my_dir_path> is the path to your working directory wherein the file 
<my_file_name> that is to be read from must exist, and which is to have a
.json extension.
"""
function fromFile(::Type{Gene}, json_stream::IOStream)::Gene
    if isopen(json_stream)
        g = JSON3.read(readline(json_stream), Gene)
    else
        msg = "The supplied JSON stream is not open."
        error(msg)
    end
    return g
end
