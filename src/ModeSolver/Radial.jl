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

#function D_coeffs(ω,a,s,m,Alm)

#end

	
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
    σ = Alm + (a*ω)^2 - 8*ω^2 + p*(2*α + γ - δ) + (1 + s - (γ + δ)/2 )*(s + (γ + δ)/2)

    D₀ = δ
    D₁ = 4*p - 2*α + γ - δ -2
    D₂ = 2*α - γ + 2
    D₃ = α*(4*p - δ) - σ
    D₄ = α*(α - γ + 1)

    return ((ζ,ξ,η),(p,α,γ,δ,σ),(D₀,D₁,D₂,D₃,D₄))
end

ParameterTransformations(P::ModeParameters) = ParameterTransformations(P.l,P.m,P.s,P.a,P.ω,P.Alm)
ParameterTransformations(l,m,n,s,a,ω,Alm,Nmax,lmax) = ParameterTransformations(l,m,s,a,ω,Alm)



function ContinuedFraction(α,β,γ,n::Int,N::Int; tol = 1e-10, Nmin = 20, rN = zero(α(1)))
    num = β(0)
    for i = 1:n
        num = β(i) - α(i-1)*γ(i)/num
    end
    lhs = num

    num = β(N)+α(N)*rN
    for i = (N):-1:(n+1)
       num = β(i-1) - α(i-1)*γ(i)/num
    end
    rhs = num - β(n)

    lhs + rhs
end


function rNCoeffs(D0,D1,D2,D3,D4)
    fp = sqrt(-D0 - D1 - D2)
    u1 = -fp
    u2 = (1/2)*(-4 - 2*D0 - D1)
    u3 = (fp)*(8 + 16*D0 + 8*(D0^2) + 12*D1 + 8*D0*D1 + (D1^2) + 8*D2 +
    4*D0*D2 - 4*D3 - 4*D4)/(8*(D0 + D1 + D2))
    u4 = (1/2)*(4 + 4*D0 + 2*D0^2 + D1 + D0*D1 - D3)
    u1,u2,u3,u4
end

function CF(l,m,n,s,a,ω,Alm; Nmax=10)
    _,_,(D₀,D₁,D₂,D₃,D₄) = ParameterTransformations(P)
    n = P.n; N = P.Nmax;
    αₙ(n) = n^2 + (D₀+1)*n + D₀
    βₙ(n) = -2*n^2 + (D₁+2)*n +D₃
    γₙ(n) = n^2 + (D₂-3)*n  + D₄ - D₂ + 2
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)
    ContinuedFraction(αₙ,βₙ,γₙ,n,N; rN = rN)
end

function CF(P)
    _,_,(D₀,D₁,D₂,D₃,D₄) = ParameterTransformations(P)
    n = P.n; N = P.Nmax;
    αₙ(n) = n^2 + (D₀+1)*n + D₀
    βₙ(n) = -2*n^2 + (D₁+2)*n +D₃
    γₙ(n) = n^2 + (D₂-3)*n  + D₄ - D₂ + 2
    u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
    rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)
    ContinuedFraction(αₙ,βₙ,γₙ,n,N; rN = rN)
end

"""
My old lentz inversion code had bugs and was slow. 
I have completely taken Leo Stein's function inside `qnm/radial.py` that solves the lentz inversion - with completely minimal edits. 
This is more stable and faster (I think even faster in julia since this is now non-allocating)
"""
function leaver_CF_inv_lentz(ω,a,s,m,Alm,n_inv; tol=1e-10, N_min=1, N_max=Inf)
	_,_,(D₀,D₁,D₂,D₃,D₄) = ParameterTransformations(0,m,s,a,ω,Alm)
    #n = 0:(n_inv+1)
	α(nᵢ::Integer) = nᵢ*nᵢ + (D₀ + 1.)*nᵢ + D₀
	β(nᵢ::Integer)  = -2.0*nᵢ*nᵢ + (D₁ + 2.)*nᵢ + D₃
	γ(nᵢ::Integer) = nᵢ*nᵢ + (D₂ - 3.)*nᵢ + D₄ - D₂ + 2.0

	#=
	nᵢ = 0:(n_inv)
    alpha = @.     nᵢ*nᵢ + (D₀ + 1.)*nᵢ + D₀
    beta  = @. -2.0*nᵢ*nᵢ + (D₁ + 2.)*nᵢ + D₃
    gamma = @.    nᵢ*nᵢ + (D₂ - 3.)*nᵢ + D₄ - D₂ + 2.

    conv1 = zero(D₀)
    for i in 1:n_inv # n_inv is not included (based on 0 indexing)
        conv1 = alpha[i] / (beta[i] - gamma[i] * conv1)
	end
	=#

	conv1 = zero(D₀)
    for i in 0:(n_inv-1)
        conv1 = α(i) / (β(i) - γ(i) * conv1)
	end
	
	tiny = 1.e-30

	# This is starting with b_0 = 0 for the infinite continued
    # fraction. I could have started with other values (e.g. b_i
    # evaluated with i=0) but then I would have had to subtract that
    # same quantity away from the final result. I don't know if this
    # affects convergence.
	
	f_old = tiny

    C_old = f_old
    D_old = zero(D₀)

	f_new = f_old;
	C_new = C_old;
	D_new = D_old;
	Delta = (C_new)*(D_new)

    conv = false

    j = 1
    n = n_inv

    while ((!conv) & (j < N_max))

        # In defining the below a, b sequences, I have cleared a fraction
        # compared to the usual way of writing the radial infinite
        # continued fraction. The point of doing this was that so both
        # terms, a(n) and b(n), tend to 1 as n goes to infinity. Further,
        # We can analytically divide through by n in the numerator and
        # denominator to make the numbers closer to 1.
        an = -(n*n + (D₀ + 1.)*n + D₀)/(n*n + (D₂ - 3.)*n + D₄ - D₂ + 2.)
        n = n + 1
        bn = (-2.0*n*n + (D₁ + 2.)*n + D₃)/(n*n + (D₂ - 3.)*n + D₄ - D₂ + 2.)

        D_new = bn + an * D_old

        if (D_new == 0)
            D_new = tiny
		end

        C_new = bn + an / C_old

        if (C_new == 0)
            C_new = tiny
		end

        D_new = 1.0/D_new
        Delta = C_new * D_new
        f_new = f_old * Delta

        if ((j > N_min) & (abs(Delta - 1.) < tol)) # converged
            conv = true
		end

		 #println(j," ",an," ", bn," ", D_new, " ", C_new, " ", f_old)

        # Set up for next iter
        j = j + 1
        D_old = D_new
        C_old = C_new
        f_old = f_new
	end

    conv2 = f_new
    ##############################
    
    return ((β(n_inv)
            - γ(n_inv) * conv1
            + γ(n_inv) * conv2), abs(Delta-1.), j-1)

	
end
