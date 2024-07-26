"""
Counter is a container for counters that appear in a struct.

An instance of type Counter is mutable; hence, its count can be changed.

mutable struct Counter
    count::Integer
end

Constructors

    c = Counter()       creates a new counter with a 'count' of 0
    c = Counter(count)  creates a new counter with a 'count' of count

Operators

    boolean: ==, ≠, <, ≤, >, ≥
    arithmetic, unary:  +, -
    arithmetic, binary: +, -, *, ÷, %, ^

Methods

    count = get(c)      returns 'count' held by Counter 'c'
    set!(c, count)      assigns 'count' to Counter 'c'
    d = copy(c)         returns a copy 'd' of Counter 'c'
    s = tostring(c)     returns a string representation 's' of Counter 'c'
"""
mutable struct Counter
    count::Integer

    # constructor

    function Counter()
        count = 0
        new(count)
    end

    function Counter(count::Integer)
        new(count)
    end
end # Counter

# boolean operators

function Base.:(==)(cL::Counter, cR::Counter)::Bool
    if cL.count == cR.count
        return true
    else
        return false
    end
end # ==

function Base.:(==)(c::Counter, count::Integer)::Bool
    if c.count == count
        return true
    else
        return false
    end
end # ==

function Base.:(==)(count::Integer, c::Counter)::Bool
    if count == c.count
        return true
    else
        return false
    end
end # ==

function Base.:≠(cL::Counter, cR::Counter)::Bool
    if cL.count ≠ cR.count
        return true
    else
        return false
    end
end # ≠

function Base.:≠(c::Counter, count::Integer)::Bool
    if c.count ≠ count
        return true
    else
        return false
    end
end # ≠

function Base.:≠(count::Integer, c::Counter)::Bool
    if count ≠ c.count
        return true
    else
        return false
    end
end # ≠

function Base.:<(cL::Counter, cR::Counter)::Bool
    if cL.count < cR.count
        return true
    else
        return false
    end
end # <

function Base.:<(c::Counter, count::Integer)::Bool
    if c.count < count
        return true
    else
        return false
    end
end # <

function Base.:<(count::Integer, c::Counter)::Bool
    if count < c.count
        return true
    else
        return false
    end
end # <

function Base.:≤(cL::Counter, cR::Counter)::Bool
    if cL.count ≤ cR.count
        return true
    else
        return false
    end
end # ≤

function Base.:≤(c::Counter, count::Integer)::Bool
    if c.count ≤ count
        return true
    else
        return false
    end
end # ≤

function Base.:≤(count::Integer, c::Counter)::Bool
    if count ≤ c.count
        return true
    else
        return false
    end
end # ≤

function Base.:>(cL::Counter, cR::Counter)::Bool
    if cL.count > cR.count
        return true
    else
        return false
    end
end # >

function Base.:>(c::Counter, count::Integer)::Bool
    if c.count > count
        return true
    else
        return false
    end
end # >

function Base.:>(count::Integer, c::Counter)::Bool
    if count > c.count
        return true
    else
        return false
    end
end # >

function Base.:≥(cL::Counter, cR::Counter)::Bool
    if cL.count ≥ cR.count
        return true
    else
        return false
    end
end # ≥

function Base.:≥(c::Counter, count::Integer)::Bool
    if c.count ≥ count
        return true
    else
        return false
    end
end # ≥

function Base.:≥(count::Integer, c::Counter)::Bool
    if count ≥ c.count
        return true
    else
        return false
    end
end # ≥

# arithmetic operators

function Base.:+(c::Counter)::Integer
    return +c.count
end # +

function Base.:-(c::Counter)::Integer
    return -c.count
end # -

function Base.:+(cL::Counter, cR::Counter)::Integer
    return cL.count + cR.count
end # +

function Base.:+(c::Counter, count::Integer)::Ingeter
    return c.count + count
end # +

function Base.:+(count::Integer, c::Counter)::Integer
    return count + c.count
end # +

function Base.:-(cL::Counter, cR::Counter)::Integer
    return cL.count - cR.count
end # -

function Base.:-(c::Counter, count::Integer)::Ingeter
    return c.count - count
end # -

function Base.:-(count::Integer, c::Counter)::Integer
    return count - c.count
end # -

function Base.:*(cL::Counter, cR::Counter)::Integer
    return cL.count * cR.count
end # *

function Base.:*(c::Counter, count::Integer)::Ingeter
    return c.count * count
end # *

function Base.:*(count::Integer, c::Counter)::Integer
    return count * c.count
end # *

function Base.:÷(cL::Counter, cR::Counter)::Integer
    return cL.count ÷ cR.count
end # ÷

function Base.:÷(c::Counter, count::Integer)::Ingeter
    return c.count ÷ count
end # ÷

function Base.:÷(count::Integer, c::Counter)::Integer
    return count ÷ c.count
end # ÷

function Base.:%(cL::Counter, cR::Counter)::Integer
    return cL.count % cR.count
end # %

function Base.:%(c::Counter, count::Integer)::Ingeter
    return c.count % count
end # %

function Base.:%(count::Integer, c::Counter)::Integer
    return count % c.count
end # %

function Base.:^(cL::Counter, cR::Counter)::Integer
    return cL.count ^ cR.count
end # ^

function Base.:^(c::Counter, count::Integer)::Ingeter
    return c.count ^ count
end # ^

function Base.:^(count::Integer, c::Counter)::Integer
    return count ^ c.count
end # ^

# methods

function Base.:(get)(c::Counter)::Integer
    return copy(c.count)
end # get

function set!(c::Counter, count::Integer)
    c.count = copy(count)
end # set!

function Base.:(copy)(c::Counter)::Counter
    count = copy(c.count)
    return Counter(count)
end # copy

function tostring(c::Counter)::String
    return String(c.count)
end # tostring
