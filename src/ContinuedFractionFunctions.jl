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
