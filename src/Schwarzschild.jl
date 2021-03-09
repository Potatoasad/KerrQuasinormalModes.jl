using CSV

function LoadLowSchwarzschildQnms()
    A = CSV.read(joinpath(@__DIR__, "SwData.csv"),NamedTuple)
    ωlist = Dict{NamedTuple,Complex{Float64}}()
    for i=1:length(A[:l])
        merge!(ωlist,Dict((l = A[:l][i], s = A[:s][i], n = A[:n][i]) => A[:ωreal][i]+im*A[:ωimag][i]))
    end
    ωlist
end

function Low_l_Guess(l,n,s)::Complex{Float64}
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

    ω = (ω₋₁*L + ω₁*(L^(-1)) + ω₃*(L^(-3))) - im*N*(ω₀ᵢ + ω₂ᵢ*(L^(-2)) - ω₄ᵢ*(L^(-4)))
    ω = (1/sqrt(27))*ω
end


function High_n_Guess(l,n,s)::Complex{Float64}
    print("Used higher n guess")
    # From H.-J. Blome, B. Mashhoon, Phys. Lett. A 110, 231 (1984)
    return (0.5 + l - im*(0.5 + n))/(3.0*sqrt(3))
end


function SchwarzschildGuess(l::Int64,n::Int64,s::Int64)::Complex{Float64}
    if l >= 3
        return High_l_Guess(l,n,s)::Complex{Float64}
    else
        return Low_l_Guess(l,n,s)::Complex{Float64}
    end
end

function SolveForParams(l,m,n,s,a,ωguess,Alm,Nmax,lmax)
    let l=l,m=m,n=n,s=s,a=a,Alm=Alm,Nmax=Nmax,lmax=lmax,ωguess=ωguess
        function f!(F,x)
            val = CF(l,m,n,s,a,x[1]+im*x[2],Alm,Nmax,lmax)
            F[1] = real(val)
            F[2] = imag(val)
        end
        sol = nlsolve(f!, [real(ωguess),imag(ωguess)],ftol = 1e-14, method = :newton)
        ωnews = sol.zero
        ωnews[1]+im*ωnews[2]
    end
end

function CF(l,m,n,s,a,ω,Alm,Nmax,lmax)
    let l=l,m=m,n=n,s=s,a=a,ω=ω,Alm=Alm,Nmax=Nmax,lmax=lmax
        _,_,(D₀,D₁,D₂,D₃,D₄) = ParameterTransformations(l,m,n,s,a,ω,Alm,Nmax,lmax)
        N = Nmax;
        αₙ(n) = n^2 + (D₀+1)*n + D₀
        βₙ(n) = -2*n^2 + (D₁+2)*n +D₃
        γₙ(n) = n^2 + (D₂-3)*n  + D₄ - D₂ + 2
        u1,u2,u3,u4 = rNCoeffs(D₀,D₁,D₂,D₃,D₄)
        rN = 1 + u1*N^(-0.5) + u2*N^(-1) + u3*N^(1.5) + u4*N^(-2)
        ContinuedFraction(αₙ,βₙ,γₙ,n,N; rN = rN)
    end
end


function SchwarzschildModes(l,n,s)
    ωguess = SchwarzschildGuess(l,n,s)::Complex{Float64}
    ωguess = Complex(ωguess)
    #println(ωguess)
    Alm = l*(l+1) - s*(s+1)
    lmax = 12; m=0; a = 0

    NmaxTrials = [100,300,500,1000,5_000,10_000,50_000,100_000]
    ω₀ = ωguess; Alm₀ = Alm
    ωnew = ωguess + 1.0;
    iterations = 1;
    while (abs(ωnew-ω₀) > 1e-10) & (iterations < 8)
        Nmax = NmaxTrials[iterations]
        ω₀ = ωnew;
        ωnew = SolveForParams(l,m,n,s,a,ωguess,Alm,Nmax,lmax)
        iterations += 1
        #println("ω diff = ",abs(ωnew-ω₀))
    end
    #println("Number of Iterations = ", iterations)
    ωnew, Alm
end
