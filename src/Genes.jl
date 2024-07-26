"""
Genes are the lowest level containers of genetic information.  Here they are
considered to contain haploid gene expressions, i.e., dominant (assigned true)
and recessive (assigned false), which are held in an object of type Expression.

A gene is a type whose datum is an expression.

struct Gene
    expression::Expression
end

Constructors

    g = Gene()                  creates a gene 'g' with random expression
    g = Gene(dominant)          creates a gene 'g' with dominant expression
    g = Gene(recessive)         creates a gene 'g' with recessive expression
    g = Gene(expression)        creates a gene 'g' with 'expression'

Operators

    ==, ≠

Methods

    e = get(g)                  returns expression 'e' held by gene 'g'
    set!(g, expression)         assigns a gene `expression` to g.expression
    c = copy(g)                 returns a copy 'c' of gene 'g'
    s = tostring(g)             returns a string representation 's' of gene 'g'
    b = isdominant(g)           returns true if g.expression == dominant
    b = isrecessive(g)          returns true if g.expression == recessive
    mutate!(g, probability)     random flip in gene expression at 'probability'
"""
struct Gene
    expression::Expression

    # constructors

    function Gene()
        gene_expression = Expression()
        new(gene_expression)
    end

    function Gene(expression::Bool)
        gene_expression = Expression(expression)
        new(gene_expression)
    end

    function Gene(expression::Expression)
        new(expression)
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

function Base.:(==)(g::Gene, e::Expression)::Bool
    if g.expression == e
        return true
    else
        return false
    end
end # ==

function Base.:(==)(e::Expression, g::Gene)::Bool
    if e == g.expression
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

function Base.:≠(g::Gene, e::Expression)::Bool
    if g.expression ≠ e
        return true
    else
        return false
    end
end # ≠

function Base.:≠(e::Expression, g::Gene)::Bool
    if e ≠ g.expression
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

function tostring(g::Gene)::String
    return tostring(g.expression)
end # tostring

function isdominant(g::Gene)::Bool
    return isdominant(g.expression)
end # isdominant

function isrecessive(g::Gene)::Bool
    return isrecessive(g.expression)
end # isrecessive

"""
mutate! addresses the occurrence of a gene's expression changing from dominant
to recessive, or vice versa, with a chance of change at probability_mutation.
"""
function mutate!(g::Gene, probability_mutation::Real)
    if probability_mutation < 0.0 || probability_mutation ≥ 1.0
        msg = "A probability of mutation must belong to unit interval [0, 1)."
        error(msg)
    end
    if probability_mutation > rand()
        if isdominant(g)
            set!(g, recessive)
        else
            set!(g, dominant)
        end
    end
    return nothing
end # mutate!
