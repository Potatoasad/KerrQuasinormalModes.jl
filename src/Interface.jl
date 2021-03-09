struct ScalarFunctions{T,F}
    coord::T
    func::F
end

### Radial Functions
struct HeunConfluentRadial
    η::Complex{Float64}; α::Complex{Float64}; ξ::Complex{Float64}; ζ::Complex{Float64}
    r₊::Float64; r₋::Float64
    coeffs::Array{Complex{Float64},1}
end


function (ψᵣ::HeunConfluentRadial)(r::Union{Float64,Complex{Float64}})
    η = ψᵣ.η;
    α = ψᵣ.α;
    ξ = ψᵣ.ξ;
    ζ = ψᵣ.ζ;
    r₊ = ψᵣ.r₊;
    r₋ = ψᵣ.r₋
    asymptoticpart = (r₊-r₋)^(α)*(im*(r-r₋))^(η-α)*(im*(r-r₊))^(ξ)*exp(ζ*r)
    x = (r-r₊)/(r-r₋)
    finalsum = 0
    for n in 1:length(ψᵣ.coeffs)
       finalsum += ψᵣ.coeffs[n]*x^(n-1)
    end
    #print(asymptoticpart)
    asymptoticpart*finalsum
end

function ComplexPlot(ψ::HeunConfluentRadial; ztopleft = 0.8 + 1.0im, zbottomright = 2.0 - 0.5im)
    img = portrait(ztopleft, zbottomright, Ψ;
        point_color = cs_d(; colormap=hsv_colors()))
    display(img)
    img
end

### Angular Functions
struct SpinWeightedSpherical
    s::Int64; l::Int64; m::Int64
end

function sterlings(n)
    return sqrt(2*pi*n)*(n/exp(1))^n
end

function SpinWeightedSphericalCalculation(z,s,l,m)
    ### Convention is that θ goes from 0 to pi, and hence cos(θ)
    ### goes from 1 to -1, which means θ/2 is in the first quadrant
    ### hence all cot terms are positive, which is why I take the
    ### positive root whereever needed.
    ### Pulled from
    ### https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics
    #Define the overall factor
    #println((l+m > 20),"  ",(l+s > 20),"  ",(l-m > 20),"  ",(l-s > 20))
    #println(((l+m > 20) | (l+s > 20) | (l-m > 20) | (l-s > 20)))
    if ((l+m > 20) | (l+s > 20) | (l-m > 20) | (l-s > 20))
        term1 = (2*l+1)*(sterlings(l+m)/sterlings(l+s))*(sterlings(l-m)/sterlings(l-s))
    else
        term1 = (2*l+1)*(factorial(l+m)/factorial(l+s))*(factorial(l-m)/factorial(l-s))
    end
    #println((l=l,m=m,s=s), (fact = factorial(l+m),), term1)
    A = sqrt(term1/(4*pi))*(-1)^m

    #The sin term right outside the summation
    sinterm = ((1-z)/2)^l

    #The summation terms
    sumterms = 0
    for r = 0:(l-s)
        consts = binomial(l-s,r)*binomial(l+s,r+s-m)*(-1)^(l-r-s)
        cotterm = ((1+z)/(1-z))^((2*r+s-m)/2)
        sumterms += consts*cotterm
    end

    return A*sinterm*sumterms
end

function (Ψ::SpinWeightedSpherical)(z)
    s = Ψ.s; l = Ψ.l; m = Ψ.m;
    if (l >= min(abs(m),abs(s))) & (l >= abs(m))
        return SpinWeightedSphericalCalculation(z,s,l,m)
    end
    return Complex(0.0)
end

function (Ψ::SpinWeightedSpherical)(z,ϕ)
    Ψ(z)*exp(im*Ψ.m*ϕ)
end

struct SpinWeightedSpheroidal
    s::Int64; l::Int64; m::Int64
    Cllʼ::Array{Complex{Float64},1}
    lmin::Int64; lmax::Int64
end

function SpinWeightedSpheroidal(s,l,m,Cllʼ)
    lmins = max(abs(s),abs(m));
    lmax = length(Cllʼ) + lmins -1
    SpinWeightedSpheroidal(s,l,m,Cllʼ,lmins, lmax)
end

function SpinWeightedSpheroidalCalculation(z,s,l,m,Cllʼ,lmin,lmax)
    N = lmax - lmin + 1
    val = Complex(0.0)
    for j = 1:N
        lʼ = j+lmin-1;
        #println(lʼ)
        val += Cllʼ[j]*SpinWeightedSphericalCalculation(z,s,lʼ,m)
    end
    val
end

function (Ψ::SpinWeightedSpheroidal)(z)
    s = Ψ.s; l = Ψ.l; m = Ψ.m;
    lmin = Ψ.lmin; lmax = Ψ.lmax; Cllʼ = Ψ.Cllʼ;
    SpinWeightedSpheroidalCalculation(z,s,l,m,Cllʼ,lmin,lmax)
end

(Ψ::SpinWeightedSpheroidal)(z,ϕ) = Ψ(z)*exp(im*Ψ.m*ϕ)


function RadialCoefficients(D₀, D₁, D₂, D₃, D₄; N = 250)
    αₙ(n) = (n+1)*(n+D₀)
    βₙ(n) = -2*n^2 + (D₁+2)*n + D₃
    γₙ(n) = (n-1)*(n+D₂-2) + D₄

    an(n) = βₙ(n)/αₙ(n) ; bn(n) = γₙ(n)/αₙ(n);

    ##Set largest rN value
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)

    ComputeSeriesFromab(an,bn; rN=rN, PreN=0)
end

function ComputeSeriesFromab(an::Function,bn::Function; N=250, rN = 0.0*im, PreN = 40)
    ##Initialize rn and fn vectors
    rₙ = zeros(Complex{Float64},N+1)
    fₙ = zeros(Complex{Float64},N+1)

    rold = rN
    ##Startup Pass for rₙ
    for n = (N+PreN):-1:(N+1)
        rnew = -bn(n)/(an(n) + rold)
        rold = rnew
    end
    rₙ[N+1] = rold;
    ##Pass for rn
    for n= N:-1:1
        rₙ[n] = -bn(n)/(an(n) + rₙ[n+1])
    end

    fₙ[1] = 1;
    ##Pass for fn
    for n = 2:(N+1)
        fₙ[n] = rₙ[n-1]*fₙ[n-1]
    end
    fₙ
end

### Combining the Mode radial and angular Information
struct QuasinormalModeFunction
    s::Int64; l::Int64; m::Int64; a::Float64
    ω::Complex{Float64}
    Alm::Complex{Float64}
    R::HeunConfluentRadial
    S::SpinWeightedSpheroidal
end


function qnmfunction(; s=-2,l=2,m=2,n=0,a=0.00)
    ω, Alm, Cllʼ = qnm(l=l,m=m,s=s,n=n, a=a)

    ((ζ,ξ,η),(p,α,γ,δ,σ),(D₀,D₁,D₂,D₃,D₄)) = ParameterTransformations(l,m,s,a,ω,Alm)
    r₊ = 1 + sqrt(1-a^2); r₋ = 1 - sqrt(1-a^2)

    ##Radial WaveFunction
    aₙ = RadialCoefficients(D₀, D₁, D₂, D₃, D₄)
    Ψᵣ = HeunConfluentRadial(η,α,ξ,ζ,r₊,r₋,aₙ)

    ##Angular WaveFunction
    Ψᵪ = SpinWeightedSpheroidal(s,l,m,Cllʼ)

    QuasinormalModeFunction(s,l,m,a,ω,Alm,Ψᵣ,Ψᵪ)
end

(Ψ::QuasinormalModeFunction)(r::Number) =  Ψ.R(r)
(Ψ::QuasinormalModeFunction)(r::Number, θ::Number) =  Ψ.R(r)*Ψ.S(cos(θ))
(Ψ::QuasinormalModeFunction)(r::Number, θ::Number, ϕ::Number) =  Ψ.R(r)*Ψ.S(cos(θ))*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(r::Number, θ::Number, ϕ::Number, t::Number) =  Ψ.R(r)*Ψ.S(cos(θ))*exp(im*Ψ.m*ϕ)*exp(-im*Ψ.ω*t)

(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r,),Tuple{Number}}) = Ψ.R(x[:r])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:θ,),Tuple{Number}}) = Ψ.S(cos(x[:θ]))
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:z,),Tuple{Number}}) = Ψ.S(x[:z])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ),Tuple{Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z),Tuple{Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ, :ϕ),Tuple{Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z, :ϕ),Tuple{Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :θ, :ϕ, :t),Tuple{Number,Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(cos(x[:θ]))*exp(im*Ψ.m*ϕ)
(Ψ::QuasinormalModeFunction)(x::NamedTuple{(:r, :z, :ϕ, :t),Tuple{Number,Number,Number,Number}}) = Ψ.R(x[:r])*Ψ.S(x[:z])*exp(im*Ψ.m*ϕ)

print("Here")
