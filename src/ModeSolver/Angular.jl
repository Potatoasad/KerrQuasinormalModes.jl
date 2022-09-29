function Fₛ(s,l,m)
if l < max(abs(s),abs(m))
 return 0.0
else
 num1 = (l+1)^2 - m^2
 num2 = (l+1)^2 - s^2
 denom = (2*l+3)*(2*l+1)*(l+1)^2
 return sqrt(num1*num2/denom)
end
end

function Gₛ(s,l,m)
 if l == zero(l)
     return 0.
 else
     num1 = l^2 - m^2
     num2 = l^2 - s^2
     denom = (4*l^2 - 1)*l^2
     return sqrt(num1*num2/denom)
 end
end

function Hₛ(s,l,m)
 if (l == zero(0)) | (s == zero(0))
     return 0.
 else
     return (-m*s)/(l*(l+1))
 end
end

Aₛ(s,l,m) = Fₛ(s,l,m)*Fₛ(s,l+1,m)
Dₛ(s,l,m) = Fₛ(s,l,m)*(Hₛ(s,l+1,m) + Hₛ(s,l,m))
Bₛ(s,l,m) = Fₛ(s,l,m)*Gₛ(s,l+1,m) + Gₛ(s,l,m)*Fₛ(s,l-1,m) + Hₛ(s,l,m)^2
Eₛ(s,l,m) = Gₛ(s,l,m)*(Hₛ(s,l-1,m) + Hₛ(s,l,m))
Cₛ(s,l,m) = Gₛ(s,l,m)*Gₛ(s,l-1,m)

function M(l,lʼ,s,m,c)
 if lʼ == l-2
     return (-c^2)*(Aₛ(s,lʼ,m))
 elseif lʼ == l-1
     return (-c^2)*(Dₛ(s,lʼ,m)) + 2*c*s*Fₛ(s,lʼ,m)
 elseif lʼ == l
     return lʼ*(lʼ + 1) - s*(s+1) - (c^2)*Bₛ(s,lʼ,m) + 2*c*s*Hₛ(s,lʼ,m)
 elseif lʼ == l+1
     return (-c^2)*(Eₛ(s,lʼ,m)) + 2*c*s*Gₛ(s,lʼ,m)
 elseif lʼ == l+2
     return (-c^2)*(Cₛ(s,lʼ,m))
 else
     return 0.0*im
 end
end

lmin(s,m) = max(abs(m),abs(s))

function ConstructM(lmax,s,m,c)
 lmins = lmin(m,s)
 N = lmax - lmins + 1
 Mₑ = zeros(typeof(c),N,N)
 #Mₑ = BandedMatrix{ComplexF64}(Zeros(N,N), (2,2))
for δl ∈ -2:2
	for l ∈ lmins:lmax
		#l = i+lmins-1; lʼ = j+lmins-1;
		i = l-lmins+1; lʼ = l + δl;
		j = lʼ-lmins+1;
		if ((i > 0) & (i ≤ N)) & ((j > 0) & (j ≤ N))
			Mₑ[i,j] = M(l,lʼ,s,m,c)
		end
	end
end
 Mₑ
end

function ComputeAₗₘCll(s,m,c, Aₗₘprev,lmax)
 Mₑ = ConstructM(lmax,s,m,c)
 F = eigen(Mₑ)
 #vals, vecs, _ = KrylovKit.eigsolve(Mₑ)
 vals = Complex.(F.values)::Array{Complex{Float64},1}
 #vals = F.values::Array{Complex{Float64},1}
 vecs = F.vectors::Array{Complex{Float64},2}
 error, j = findmin(map(Anew -> abs2(Anew - Aₗₘprev), vals))
 vals[j], vecs[:,j], error
 #vals[j], vecs[j], error
end

function ComputeAₗₘ(s,m,c, Aₗₘprev,lmax)
 Mₑ = ConstructM(lmax,s,m,c)
 vals = Complex.(eigvals(Mₑ))::Array{Complex{Float64},1}
 #vals, _, _ = KrylovKit.eigsolve(Mₑ)
 #vals = Complex.(F.values)::Array{Complex{Float64},1}
 #vals = F.values::Array{Complex{Float64},1}
 #vecs = F.vectors::Array{Complex{Float64},2}
 error, j = findmin(map(Anew -> abs2(Anew - Aₗₘprev), vals))
 vals[j]
end

	
function schw_Alm(s,l,m)
	return l*(l+1) - s*(s+1)
end
