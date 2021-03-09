###Define a structure to hold all our parameters

mutable struct ModeParameters
    l::Int64
    m::Int64
    n::Int64
    s::Int64
    a::Float64
    ω::Complex{Float64}
    Alm::Complex{Float64}
    Nmax::Int64
    lmax::Int64
end

ModeParameters(NT::NamedTuple) = ModeParameters(NT[:l],NT[:m],NT[:n],NT[:s],NT[:a],NT[:ω],NT[:Alm],NT[:Nmax],NT[:lmax])

##Parameter Transformation
function ParameterTransformations(l,m,s,a,ω,Alm)
    M = 1
    r₊ = M + sqrt(M^2 - a^2)
    r₋ = M - sqrt(M^2 - a^2)

    σ₊ = (2*ω*M*r₊ - m*a)/(r₊-r₋)
    σ₋ = (2*ω*M*r₋ - m*a)/(r₊-r₋)

    ζ₊ = im*ω
    ξ₋ = (-s-(s+2*im*σ₊))/2
    η₋ = -s+im*σ₋

    ζ = ζ₊
    ξ = ξ₋
    η = η₋

    p = (r₊ - r₋)*(ζ/2)
    α = 1 + s + ξ + η - 2*ζ + s*im*ω/ζ
    γ = 1 + s + 2*η
    δ = 1 + s + 2*ξ
    σ = Alm + (a*ω)^2 - 8*ω^2 + p*(2*α + γ - δ) + (1 + s - (γ + δ)/2 )*(s + (γ + δ)/2 )

    D₀ = δ
    D₁ = 4*p - 2*α + γ - δ -2
    D₂ = 2*α - γ + 2
    D₃ = α*(4*p - δ) - σ
    D₄ = α*(α - γ + 1)

    return ((ζ,ξ,η),(p,α,γ,δ,σ),(D₀,D₁,D₂,D₃,D₄))
end

ParameterTransformations(P::ModeParameters) = ParameterTransformations(P.l,P.m,P.s,P.a,P.ω,P.Alm)
ParameterTransformations(l,m,n,s,a,ω,Alm,Nmax,lmax) = ParameterTransformations(l,m,s,a,ω,Alm)
