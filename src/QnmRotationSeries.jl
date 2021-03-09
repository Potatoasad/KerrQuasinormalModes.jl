struct LinearInterpolant{T,N}
    y::Array{T,N}
    Δx::Float64
end


function (f::LinearInterpolant)(a::Float64)
    Δx = f.Δx
    s = a/(Δx);
    n₁ = floor(Int64,s);
    n₂ = ceil(Int64,s);
    Δn = (s - n₁)
    (1-Δn)*f.y[n₁+1] + Δn*f.y[n₂+1]
end

function qnms(;l=0,m=0,n=0,s=0,amax=0.99, ϵ = 0.01)
    as, ωs, Almss, Cllss, P, ϵ = GetModes(l,m,n,s;amax = amax,ϵ = ϵ)
    ωfunc = LinearInterpolant(ωs,ϵ)
    Afunc = LinearInterpolant(Almss,ϵ)
    Cllss = LinearInterpolant(Cllss,ϵ)
    ωfunc, Afunc, Cllss
end


function qnm(;l=0,m=0,n=0,s=0,a=0.0, ϵ = 0.01)
     _, ωs, Almss, Cllss, P, _ = GetModes(l,m,n,s;amax = a,ϵ = ϵ)
    #ωs[end], Almss[end], Cllss[:,end], P
    Alms, Cll = ComputeAₗₘ(P.s,P.m,P.a*P.ω, P.Alm,P.lmax)
    ωs[end], Alms, Cll, P
end
