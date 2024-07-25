# Admissible gene expressions include:

const dominant  = true    # Gene expression whenever a gene is dominant.
const recessive = false   # Gene expression whenever a gene is recessive.

"""
Expression is a container for gene expression.  Haploid genes are considered,
i.e., they can be either dominant (assigned true) or recessive (assigned false).

An instance of type Expression is mutable; hence, its expression can be changed.

mutable struct Expression
    expression::Bool
end

Constructors

    e = Expression()            creates a  random   expression 'e'
    e = Expression(dominant)    creates a dominant  expression 'e'
    e = Expression(recessive)   creates a recessive expression 'e'

Operators

    ==, ≠

Methods

    expression = get(e)     returns 'expression' held by Expression 'e'
    set!(e, expression)     assigns 'expression' to Expression 'e'
    c = copy(e)             returns a copy 'c' of Expression 'e'
    s = tostring(e)         returns string representation 's' of Expression 'e'
    b = isdominant(e)       returns true if e.expression == dominant
    b = isrecessive(e)      returns true if e.expression == recessive
"""
mutable struct Expression
    expression::Bool

    # constructor

    function Expression()
        value = rand()    # Returns a random real from the interval [0, 1).
        if value < 0.5
            expression = recessive
        else
            expression = dominant
        end
        new(expression)
    end
end # Expression

# operators

function Base.:(==)(eL::Expression, eR::Expression)::Bool
    if eL.expression == eR.expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(e::Expression, expression::Bool)::Bool
    if e.expression == expression
        return true
    else
        return false
    end
end # ==

function Base.:(==)(expression::Bool, e::Expression)::Bool
    if expression == e.expression
        return true
    else
        return false
    end
end # ==

function Base.:≠(eL::Expression, eR::Expression)::Bool
    if eL.expression ≠ eR.expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(e::Expression, expression::Bool)::Bool
    if e.expression ≠ expression
        return true
    else
        return false
    end
end # ≠

function Base.:≠(expression::Bool, e::Expression)::Bool
    if expression ≠ e.expression
        return true
    else
        return false
    end
end # ≠

# methods

function Base.:(get)(e::Expression)::Bool
    return copy(e.expression)
end # get

function set!(e::Expression, expression::Bool)
    e.expression = copy(expression)
    return nothing
end # set!

function Base.:(copy)(e::Expression)::Expression
    expression = copy(e.expression)
    return Expression(expression)
end # copy

function tostring(e::Expression)::String
    if e.expression == dominant
        str = "1"
    else
        str = "0"
    end
    return str
end # tostring

function isdominant(e::Expression)::Bool
    if e.expression == dominant
        return true
    else
        return false
    end
end # isdominant

function isrecessive(e::Expression)::Bool
    if e.expression == recessive
        return true
    else
        return false
    end
end # isrecessive
