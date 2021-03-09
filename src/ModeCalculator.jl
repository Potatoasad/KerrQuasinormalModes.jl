function GetModes(l,m,n,s; amax = 0.99, ϵ = 0.01, Nmax = 300)
    ##Initialize
    as = [a for a in 0.0:ϵ:amax]
    ωs = zero(im*(as))
    Almss = zero(im*(as))
    Cllss = zero([zero(22) for i in 1:length(as)])
    NCsize = 22; Nas = 100
    Cllss = zeros(Complex{Float64},NCsize,Nas)

    #Get Schwarzschild Modes
    ωsch, Almsch = SchwarzschildModes(l,n,s)

    ϵs = 0.0000001 #Finite difference step size

    P = ModeParameters((l=l,s=s,m=m,n=n,a=0.0,ω = ωsch,Alm = Almsch, Nmax = Nmax, lmax = 10))

    NmaxTrials = [100,300,500,1000,5_000,10_000,50_000]
    #NmaxTrials = [10,15,30,50,100,500,1000,5000]
    lmaxTrials = [10,12,15,15,20,22]

    for i =1:length(as)
        P.a = as[i]
        ω₀ = P.ω; Alm₀ = P.Alm
        ωnew = P.ω + 1.0; Almnew = Alm₀ + 1.0;
        iterations = 1;
        while ((abs(ωnew-ω₀) > 1e-10) | (abs(Almnew-Alm₀) > 1e-10)) & (iterations < 6)
            P.Nmax = NmaxTrials[iterations]
            P.lmax = lmaxTrials[iterations]
            ω₀ = P.ω; Alm₀ = P.Alm;
            Newton!(P,ϵs)
            ωnew = P.ω; Almnew = P.Alm
            iterations += 1
        end

        ωs[i] = P.ω
        Almss[i] = P.Alm
        _,Cll,_ = ComputeAₗₘ(P.s,P.m,P.a*P.ω, P.Alm,P.lmax)
        for j in length(Cll)
            for i in Nas
                Cllss[j,i]  = Cll[j]
            end
        end


    end
    as,ωs,Almss,Cllss,P, ϵ
end
