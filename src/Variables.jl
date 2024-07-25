"""
Variable is a container for variables that appear in a struct.

An instance of type Variable is mutable; hence, its value can be changed.

mutable struct Variable
    value::Real
end

Constructors

    v = Variable()      creates a new variable with a 'value' of 0.0
    v = Variable(value) creates a new variable with a 'value' of value

Operators

    boolean: ==, ≠, <, ≤, >, ≥
    arithmetic, unary:  +, -
    arithmetic, binary: +, -, *, /, ^

Methods

    value = get(v)      returns 'value' held by Variable 'v'
    set!(v, value)      assigns 'value' to Variable 'v'
    c = copy(v)         returns a copy 'c' of Variable 'v'
    s = tostring(v)     returns a string representation 's' of Variable 'v'
"""
mutable struct Variable
    value::Real

    # constructor

    function Variable()
        value = 0.0
        new(value)
    end
end # Variable

# boolean operators

function Base.:(==)(vL::Variable, vR::Variable)::Bool
    if vL.value == vR.value
        return true
    else
        return false
    end
end # ==

function Base.:(==)(v::Variable, value::Real)::Bool
    if v.value == value
        return true
    else
        return false
    end
end # ==

function Base.:(==)(value::Real, v::Variable)::Bool
    if value == v.value
        return true
    else
        return false
    end
end # ==

function Base.:≠(vL::Variable, vR::Variable)::Bool
    if vL.value ≠ vR.value
        return true
    else
        return false
    end
end # ≠

function Base.:≠(v::Variable, value::Real)::Bool
    if v.value ≠ value
        return true
    else
        return false
    end
end # ≠

function Base.:≠(value::Real, v::Variable)::Bool
    if value ≠ v.value
        return true
    else
        return false
    end
end # ≠

function Base.:<(vL::Variable, vR::Variable)::Bool
    if vL.value < vR.value
        return true
    else
        return false
    end
end # <

function Base.:<(v::Variable, value::Real)::Bool
    if v.value < value
        return true
    else
        return false
    end
end # <

function Base.:<(value::Real, v::Variable)::Bool
    if value < v.value
        return true
    else
        return false
    end
end # <

function Base.:≤(vL::Variable, vR::Variable)::Bool
    if vL.value ≤ vR.value
        return true
    else
        return false
    end
end # ≤

function Base.:≤(v::Variable, value::Real)::Bool
    if v.value ≤ value
        return true
    else
        return false
    end
end # ≤

function Base.:≤(value::Real, v::Variable)::Bool
    if value ≤ v.value
        return true
    else
        return false
    end
end # ≤

function Base.:>(vL::Variable, vR::Variable)::Bool
    if vL.value > vR.value
        return true
    else
        return false
    end
end # >

function Base.:>(v::Variable, value::Real)::Bool
    if v.value > value
        return true
    else
        return false
    end
end # >

function Base.:>(value::Real, v::Variable)::Bool
    if value > v.value
        return true
    else
        return false
    end
end # >

function Base.:≥(vL::Variable, vR::Variable)::Bool
    if vL.value ≥ vR.value
        return true
    else
        return false
    end
end # ≥

function Base.:≥(v::Variable, value::Real)::Bool
    if v.value ≥ value
        return true
    else
        return false
    end
end # ≥

function Base.:≥(value::Real, v::Variable)::Bool
    if value ≥ v.value
        return true
    else
        return false
    end
end # ≥

# arithmetic operators

function Base.:+(v::Variable)::Real
    return +v.value
end # +

function Base.:-(v::Variable)::Real
    return -v.value
end # -

function Base.:+(vL::Variable, vR::Variable)::Real
    return vL.value + vR.value
end # +

function Base.:+(v::Variable, value::Real)::Real
    return v.value + value
end # +

function Base.:+(value::Real, v::Variable)::Real
    return value + v.value
end # +

function Base.:-(vL::Variable, vR::Variable)::Real
    return vL.value - vR.value
end # -

function Base.:-(v::Variable, value::Real)::Real
    return v.value - value
end # -

function Base.:-(value::Real, v::Variable)::Real
    return value - v.value
end # -

function Base.:*(vL::Variable, vR::Variable)::Real
    return vL.value * vR.value
end # *

function Base.:*(v::Variable, value::Real)::Real
    return v.value * value
end # *

function Base.:*(value::Real, v::Variable)::Real
    return value * v.value
end # *

function Base.:/(vL::Variable, vR::Variable)::Real
    return vL.value / vR.value
end # /

function Base.:/(v::Variable, value::Real)::Real
    return v.value / value
end # /

function Base.:/(value::Real, v::Variable)::Real
    return value / v.value
end # /

function Base.:^(vL::Variable, vR::Variable)::Real
    return vL.value ^ vR.value
end # ^

function Base.:^(v::Variable, value::Real)::Real
    return v.value ^ value
end # ^

function Base.:^(value::Real, v::Variable)::Real
    return value ^ v.value
end # ^

# methods

function Base.:(get)(v::Variable)::Real
    return copy(v.value)
end # get

function set!(v::Variable, value::Real)
    v.value = copy(value)
end # set!

function Base.:(copy)(v::Varaible)::Variable
    value = copy(v.value)
    return Variable(value)
end # copy

function Base.:(deepcopy)(v::Variable)::Variable
    value = deepcopy(v.value)
    return Variable(value)
end # deepcopy

function tostring(v::Variable; significant_figures::Integer)::String
    if significant_figures ≤ 2
        s = @sprintf "%.1e" v.value;
    elseif significant_figures == 3
        s = @sprintf "%.2e" v.value;
    elseif significant_figures == 4
        s = @sprintf "%.3e" v.value;
    elseif significant_figures == 5
        s = @sprintf "%.4e" v.value;
    elseif significant_figures == 6
        s = @sprintf "%.5e" v.value;
    else
        s = @sprintf "%.6e" v.value;
    end
    return s
end # tostring
