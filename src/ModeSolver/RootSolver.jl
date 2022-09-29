@with_kw struct Mode{K <: Integer}
	s::K            = -2
	l::K            = 2
	m ::K           = 2
	n::K            = 2
end

@with_kw struct RootOptions{T <: Real, K <: Integer}
	mode::Mode{K}     = Mode(s=-2,l=2,m=2,n=2)
	l_max::K         = 20
	tol::T          = √(eps(Float64))
	cf_tol::T       = 1e-10
	N_min::K        = 1
	N_max::K        = 300
end

# Root polish given a guess
function (R::RootOptions)(a,ω_guess,Alm_guess)
	@unpack mode,l_max,tol,cf_tol,N_min,N_max = R
	@unpack s,m,l,n = mode
	F! = function solver!(Xnew,x)
			ω = x[1]+im*x[2]
			Alm = ComputeAₗₘ(s,m,a*ω, Alm_guess,l_max)
			ω,_,_ = leaver_CF_inv_lentz(ω,a,s,m,Alm,n; tol=cf_tol, N_min=N_min, N_max=N_max)
			Xnew[1] = real(ω);
			Xnew[2] = imag(ω);
		end	
	result = nlsolve(F!, [real(ω_guess),imag(ω_guess)])
	ω_new = result.zero[1] + result.zero[2]*(1.0*im)
	Alm_new = ComputeAₗₘ(s,m,a*ω_new, Alm_guess,l_max)
	ω_new,Alm_new,result.zero
end
