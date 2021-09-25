#### Creating a small symbolic system for linear sums of the qnmfunctions
#### This is so we can write the derivative as a sum over multiple qnmfunctions

struct LinearCombinationOf{T <: CallableAtom} <: Callable
    dict::Dict{T,Complex{Float64}}
    constant::Complex{Float64}
end

LinearCombinationOf(x;c = Complex(0.0)) = LinearCombinationOf(x,c)

#(x::IdentityFunction)(s) = Complex(1.0)
#(x::HeunConfluentRadial)(s) = s^2 + 2s
(x::LinearCombinationOf{T})(s...) where T = sum(k(s...)*v for (k,v) ∈ x.dict) + x.constant

function postprocess(x::LinearCombinationOf{T}) where T
    if all(x == zero(x) for x ∈ values(x.dict))
        return first(x.dict)
    end
    LinearCombinationOf(filter(a->a.second!=zero(a.second),x.dict),x.constant)
end

import Base.+, Base.-, Base.*
#const HCRNum = Union{HeunConfluentRadial, Number}
const LC = LinearCombinationOf{HeunConfluentRadial}
const HCR = HeunConfluentRadial
const LCHCR = Union{HeunConfluentRadial,LinearCombinationOf{HeunConfluentRadial}}

GetDict(x::Number) = Dict()
GetDict(x::CallableAtom; val = 1.0) = Dict(x => convert(Complex{Float64},val))
GetDict(x::LinearCombinationOf{T}) where T = x.dict

constant(x::Number) = x
constant(x::CallableAtom) = Complex(0.0)
constant(x::LinearCombinationOf{T}) where T = x.constant


## Overload Plus

function (+)(x::CallableAtom)
    LinearCombinationOf(GetDict(x))
end

function (+)(x::Callable,y::Callable)
    z = LinearCombinationOf(merge(+,GetDict(x),GetDict(y)),constant(x)+constant(y))
    postprocess(z)
end

function (+)(x::Number,y::Callable)
    LinearCombinationOf(GetDict(y),constant(x)+constant(y))
end

(+)(y::Callable, x::Number) = x+y

## Overload Minus

function (-)(x::CallableAtom)
    LinearCombinationOf(GetDict(x; val = -1.0))
end

function (-)(x::Callable,y::Callable)
    z = LinearCombinationOf(merge(+,GetDict(x),Dict((k,-v) for (k,v) ∈ GetDict(y))),constant(x)+constant(y))
    postprocess(z)
end

function (-)(x::Number,y::Callable)
    LinearCombinationOf(GetDict(y),constant(x)-constant(y))
end

function (-)(y::Callable,x::Number)
    LinearCombinationOf(GetDict(y),constant(y)-constant(x))
end

## Overload Times

function (*)(x::Number,y::Callable)
    ydict = GetDict(y);
    z = LinearCombinationOf(Dict((k,v*x) for (k,v) ∈ GetDict(y)),constant(y)*x)
    postprocess(z)
end

(*)(y::Callable,x::Number) = x*y

## Overload Show

function Base.show(io::IO, ψ::LinearCombinationOf{T}) where {T}
    N = length(ψ.dict); n = 1;
    print(io,"   ")
    for (k,v) in ψ.dict
        print(io,"(",v,")")
        print(io,k)
        if n < N
            println(io," + ")
        end
        n += 1
    end
end
