struct SpinSequence{L <: AbstractArray, T <: SplineInterpolation,K<:Integer,L1<:Number,T1 <: Real, K1 <: Integer}
	mode::Mode{K}
	a::L
	ωᵣ::T
	ωᵢ::T
	Almᵣ::T
	Almᵢ::T
	sequence_options::SpinSequenceOptions{L1}
	root_options::RootOptions{T1, K1}
end

function (x::SpinSequence)(a_val)
	ω_new = x.ωᵣ(a_val) + im*x.ωᵢ(a_val)
	Alm_new = x.Almᵣ(a_val) + im*x.Almᵢ(a_val)
	Alm, Cll = ComputeAₗₘCll(x.mode.s,x.mode.m,a_val*ω_new, Alm_new,x.root_options.l_max)
	ω_new, Alm, Cll
end

function SpinSequence(s,l,m,n; a_max=0.99)
	ω0, Alm0 = SchwarszchildValues(s,l,n)
	mode = Mode(s=s,l=l,m=m,n=n)
	sequence_options = SpinSequenceOptions(ω_0=ω0, Alm_0=Alm0,a_max=a_max)
	root_options = RootOptions(mode=mode)
	a_s, ωs, Alms = MakeSpinSequence(
		sequence_options,
		root_options
	)
	ωᵣ = interpolate(a_s,real.(ωs), BSplineOrder(3))
	ωᵢ = interpolate(a_s,imag.(ωs), BSplineOrder(3))
	Almᵣ = interpolate(a_s,real.(Alms), BSplineOrder(3))
	Almᵢ = interpolate(a_s,imag.(Alms), BSplineOrder(3))
	SpinSequence(mode, a_s, ωᵣ, ωᵢ, Almᵣ, Almᵢ, sequence_options, root_options)
end

SpinSequence(m::Mode; a_max=0.99) = SpinSequence(m.s,m.l,m.m,m.n; a_max=a_max)

function SpinSequence(y::RootOptions, x::SpinSequenceOptions)
	mode = y.mode
	@unpack s,l,m,n = mode
	a_s, ωs, Alms = MakeSpinSequence(x,y)
	ωᵣ = interpolate(a_s,real.(ωs), BSplineOrder(3))
	ωᵢ = interpolate(a_s,imag.(ωs), BSplineOrder(3))
	Almᵣ = interpolate(a_s,real.(Alms), BSplineOrder(3))
	Almᵢ = interpolate(a_s,imag.(Alms), BSplineOrder(3))
	SpinSequence(mode, a_s, ωᵣ, ωᵢ, Almᵣ, Almᵢ, sequence_options, root_options)
end
