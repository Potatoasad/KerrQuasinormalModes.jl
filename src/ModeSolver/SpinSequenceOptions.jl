@with_kw struct SpinSequenceOptions{T<:Number}
	a_0::Float64 = 0.0
	a_max::Float64 = 0.9995
	delta_a_min::Float64 = 1e-5
	delta_a_max::Float64 = 5e-3
	Alm_0::T
	ω_0::T
end

function SchwarszchildValues(s,l,n)
	#ω = 0.37367168441804177	-0.08896231568893546*im
	#Alm = l*(l+1) - s*(s+1) + 0.0*im
	#ω,Alm
	SchwarzschildQNM(s,l,n)
end

function ∂²_edge(xdata,ydata)
	Y = interpolate(xdata, ydata, BSplineOrder(3))
	∂²Y = diff(diff(Y))
	Y(xdata[end]), ∂²Y(xdata[end])
end

function next_δa_and_predictions(i,a,ωᵣ,ωᵢ,Almᵣ,Almᵢ,DELTA_MAX, DELTA_MIN)
	last3 = (i-3):(i-1)
	ωᵣ_new, ∂²ωᵣ = ∂²_edge(a[last3],ωᵣ[last3])
	ωᵢ_new, ∂²ωᵢ = ∂²_edge(a[last3],ωᵢ[last3])
	Almᵣ_new, ∂²Almᵣ = ∂²_edge(a[last3],Almᵣ[last3])
	Almᵢ_new, ∂²Almᵢ = ∂²_edge(a[last3],Almᵢ[last3])
	denom = max(∂²ωᵣ^2 + ∂²ωᵢ^2, ∂²Almᵣ^2 + ∂²Almᵢ^2)
	δa = 0.05/√(√(abs(denom)))
	δa = max(δa, DELTA_MIN)
	δa = min(δa, DELTA_MAX)
	return δa, ωᵣ_new + im*ωᵢ_new, Almᵣ_new + im*Almᵢ_new
end

function MakeSpinSequence(x::SpinSequenceOptions, R::RootOptions)
	@unpack a_0, delta_a_min, delta_a_max, Alm_0, ω_0 = x

	#Preallocate an array that's about 10 times that from the min step
	N = 3
	#a = LinRange(x.a_0, x.a_max, N) |> collect;
	a = [x.a_0, x.a_0 + x.delta_a_max, x.a_0 + 2*x.delta_a_max]
	ωᵣ = zeros(N); ωᵢ = zeros(N);
	Almᵣ = zeros(N); Almᵢ = zeros(N);

	i = 1; #iteration number
	# Set standard to Schwarszchild Value
	ωᵣ[i] = real(ω_0); ωᵢ[i] = imag(ω_0);
	Almᵣ[i] = real(Alm_0); Almᵢ[i] = imag(Alm_0);

	while (i < 3)
		i += 1
		ω,Alm,_ = R(a[i],ωᵣ[i-1] + im*ωᵢ[i-1],Almᵣ[i-1] + im*Almᵢ[i-1])
		ωᵣ[i] = real(ω); ωᵢ[i] = imag(ω);
		Almᵣ[i] = real(Alm); Almᵢ[i] = imag(Alm);
	end
	
	while (a[i] < x.a_max)
		i += 1
		# Choose the right next step
		δa, ω_guess, Alm_guess = next_δa_and_predictions(i,a,ωᵣ,ωᵢ,Almᵣ,Almᵢ,x.delta_a_max, x.delta_a_min)
		if (a[i-1] + δa) > x.a_max
			break
		end
		push!(a,a[i-1] + δa)

		# Compute the new thing
		ω,Alm,_ = R(a[i],ω_guess, Alm_guess)
		push!(ωᵣ, real(ω));  push!(ωᵢ, imag(ω));
		push!(Almᵣ, real(Alm)); push!(Almᵢ, imag(Alm));
	end

	if a[end] != x.a_max
		push!(a,x.a_max)
		δa, ω_guess, Alm_guess = next_δa_and_predictions(i,a,ωᵣ,ωᵢ,Almᵣ,Almᵢ,x.delta_a_max, x.delta_a_min)

		# Compute the new thing
		ω,Alm,_ = R(a[i],ω_guess, Alm_guess)
		push!(ωᵣ, real(ω));  push!(ωᵢ, imag(ω));
		push!(Almᵣ, real(Alm)); push!(Almᵢ, imag(Alm));
	end

	a, ωᵣ .+ im.*ωᵢ, Almᵣ .+ im.*Almᵢ 
end
