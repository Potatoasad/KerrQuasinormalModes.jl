function LoadLowSchwarzschildQnms()
    A = CSV.read(joinpath(@__DIR__, "SwData.csv"),NamedTuple)
    ωlist = Dict{NamedTuple,Complex{Float64}}()
    for i=1:length(A[:l])
        merge!(ωlist,Dict((l = A[:l][i], s = A[:s][i], n = A[:n][i]) => A[:ωreal][i]+im*A[:ωimag][i]))
    end
    ωlist
end

function Lookup_Guess(l,n,s)::Complex{Float64}
    SchwMODESLOWWℓ[(l=l,s=s,n=n)]
end

SchwMODESLOWWℓ = LoadLowSchwarzschildQnms()

function High_l_Guess(l,n,s)
    #Taken Directly from Eq: 17-22 of https://arxiv.org/pdf/0908.0329.pdf
    L = l+0.5
    N = n+0.5
    β = 1-s^2

    ω₋₁ = 1
    ω₀ᵢ = 1
    ω₁ = β/3 - (5/36)*N^2 - (115/432)
    ω₂ᵢ = β/9 + (235/3888)*N^2 - (1415/15552)
    ω₃ = -(1/27)*β^2 + β*(204*N^2 +211)/3888 + (854160*N^4 - 1664760*N^2 - 776939)/(40310784)
    ω₄ᵢ = (1/27)*β^2 + (1100*N^2 - 2719)/(46656) + (11273136*N^4 - 52753800*N^2 + 66480535)/(2902376448)

#    ω = (ω₋₁*L + ω₁*(L^(-1)) + ω₃*(L^(-3))) - im*N*(ω₀ᵢ + ω₂ᵢ*(L^(-2)) - ω₄ᵢ*(L^(-4)))
	    ω = (ω₋₁*L + ω₁*(L^(-1)) + ω₃*(L^(-3))) - im*N*(ω₀ᵢ + ω₂ᵢ*(L^(-2)) + ω₄ᵢ*(L^(-4)))
    ω = (1/sqrt(27))*ω
end

#=
function High_n_Guess(l,n,s)::Complex{Float64}
    #print("Used higher n guess")
    # From H.-J. Blome, B. Mashhoon, Phys. Lett. A 110, 231 (1984)
    return (0.5 + l - im*(0.5 + n))/(3.0*sqrt(3))
end
=#
function High_n_Guess(l,n,s)::Complex{Float64}
    k = log(3.0)/(8*π)
    kappa = 0.25
    return k - 1.0*im * kappa * (n + 0.5)
end

function SchwarzschildGuess(s::Int64, l::Int64,n::Int64)::Complex{Float64}
	if (abs(s) ≤ l) & (l ≤ 2) & (n ≤ 10)
		return Lookup_Guess(l,n,s)::Complex{Float64}
	end
	if (( n > 3 ) & (n >= 2*l))
        return High_n_Guess(l,n,s)
    else
        return High_l_Guess(l,n,s)
	end
end

function SchwarzschildQNM_Cll(s,l,n)
	RO1 = RootOptions(mode=Mode(s=s,l=l,n=n))
	ω_new, Alm_new = RO1(0.0, SchwarzschildGuess(s,l,n), l*(l+1) - s*(s+1))
	Alm_new, Cll = ComputeAₗₘCll(s,0,0, Complex(Alm_new),RO1.l_max)
	ω_new, Alm_new, Cll
end

function SchwarzschildQNM(s,l,n)
	RO1 = RootOptions(mode=Mode(s=s,l=l,n=n))
	ω_new, Alm_new = RO1(0.0, SchwarzschildGuess(s,l,n), l*(l+1) - s*(s+1))
	Alm_new_2, Cll = ComputeAₗₘCll(s,0,Complex(0.0), Complex(Alm_new),RO1.l_max)
	ω_new, Alm_new_2
end