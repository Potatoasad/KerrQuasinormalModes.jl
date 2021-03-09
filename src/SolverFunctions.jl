
function FiniteDifferenceAtParameters(P,ϵs)
    y = CF(P)
    ωᵢ = P.ω
    Almᵢ = P.Alm

    P.ω = ωᵢ + ϵs+im*0
    P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), Almᵢ,P.lmax)
    yʼr = (CF(P) - y)/ϵs

    P.ω = ωᵢ + im*ϵs;
    P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), Almᵢ,P.lmax)
    yʼi = (CF(P) - y)/ϵs

    P.ω = ωᵢ
    P.Alm = Almᵢ
    return yʼr,yʼi,y
end

function Newton!(P,ϵs)
    ω₀ = P.ω
    Alm₀ = P.Alm
    function f!(F,x)
        ωᵢ = P.ω
        Almᵢ = P.Alm
        P.ω = x[1]+im*x[2]
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
        val = CF(P)
        P.ω = ωᵢ
        P.Alm = Almᵢ
        F[1] = real(val)
        F[2] = imag(val)
    end

    function j!(J,x)
        ω = x[1]+im*x[2]
        P.ω = ω
        ydr,ydi,y = FiniteDifferenceAtParameters(P,ϵs)
        J[1,1] = real(ydr)
        J[1,2] = real(ydi)
        J[2,1] = imag(ydr)
        J[2,2] = imag(ydi)
    end

    function fj!(F,J,x)
        ωᵢ = P.ω
        Almᵢ = P.Alm
        P.ω = x[1]+im*x[2]
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), Almᵢ,P.lmax)
        val = CF(P)
        F[1] = real(val)
        F[2] = imag(val)

        P.ω = x[1]+im*x[2] + ϵs+im*0
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), Almᵢ,P.lmax)
        yʼr = (CF(P) - val)/ϵs

        P.ω = x[1]+im*x[2] + im*ϵs
        P.Alm,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), Almᵢ,P.lmax)
        yʼi = (CF(P) - val)/ϵs

        P.ω = ωᵢ
        P.Alm = Almᵢ

        J[1,1] = real(yʼr )
        J[1,2] = real(yʼi)
        J[2,1] = imag(yʼr )
        J[2,2] = imag(yʼi)
    end

    x0 = [real(ω₀),imag(ω₀)]
    df = OnceDifferentiable(f!, j!, fj!,x0,similar(x0))
    sol = nlsolve(df, x0, ftol = 1e-14)

    ωnew = sol.zero
    P.ω = ωnew[1] + im*ωnew[2]
    Almnew,_,_ = ComputeAₗₘ(P.s,P.m,P.a*(P.ω), P.Alm,P.lmax)
    P.Alm = Almnew
    P
end
