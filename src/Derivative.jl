## Define the derivative of a qnm function

#Define the derivative of a HeunConfluentRadial function
function ∂ᵣ(Ψᵣ::HeunConfluentRadial)
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊
    r₋ = ψᵣ.r₋
    aₙ = Ψᵣ.coeffs
    """Add a sum over different copies of qnm with some
    change"""
    Ψη = HeunConfluentRadial(η-1,α,ξ,ζ,r₊,r₋,aₙ)
    Ψξ = HeunConfluentRadial(η,α,ξ-1,ζ,r₊,r₋,aₙ)
    aₙshift = circshift(aₙ,-1)
    aₙshift[end] = 0
    nn = 1:length(aₙshift)
    aₙshift = aₙshift .*nn
    Ψaₙ = HeunConfluentRadial(η-1,α+1,ξ,ζ,r₊,r₋,aₙshift)
    (im*(η-α))*Ψη + (im*ξ)*Ψξ + ζ*ψᵣ - Ψaₙ
end

#Define the derivation of the qnm function
function ∂ᵣ(Ψ::QuasinormalModeFunction)
    s = Ψ.s; l = Ψ.l; m = Ψ.m; n = Ψ.n; a = Ψ.a; ω = Ψ.ω; Alm = Ψ.Alm;
    ψᵣ = Ψ.R
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊
    r₋ = ψᵣ.r₋
    aₙ = ψᵣ.coeffs
    """Add a sum over different copies of qnm with some
    change"""
    Ψη = HeunConfluentRadial(η-1,α,ξ,ζ,r₊,r₋,aₙ)
    Ψηf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψη,Ψ.S)
    Ψξ = HeunConfluentRadial(η,α,ξ-1,ζ,r₊,r₋,aₙ)
    Ψξf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψξ,Ψ.S)
    aₙshift = circshift(aₙ,-1)
    aₙshift[end] = 0.0*im
    nn = 1:length(aₙshift)
    #print("here")
    aₙshift = aₙshift .*nn
    Ψaₙ = HeunConfluentRadial(η-1,α+1,ξ,ζ,r₊,r₋,aₙshift)
    Ψaₙf = QuasinormalModeFunction(s,l,m,n,a,ω,Alm,Ψaₙ,Ψ.S)
    (im*(η-α))*Ψηf + (im*ξ)*Ψξf + ζ*Ψ - Ψaₙf
end


function ∂ᵣ(Ψ::LinearCombinationOf{T}) where T
    sum(v*∂ᵣ(k) for (k,v) in Ψ.dict)
end

export ∂ᵣ
