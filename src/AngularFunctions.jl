
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

function M(l,lʼ,s,m,c)::Complex{Float64}
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
 Mₑ = zeros(Complex{Float64},N,N)
 for i = 1:N, j = 1:N
     l = i+lmins-1; lʼ = j+lmins-1;
     if abs(l-lʼ) <= 2
         Mₑ[i,j] = M(l,lʼ,s,m,c)
     end
 end
 Mₑ
end

function ComputeAₗₘ(s,m,c, Aₗₘprev,lmax)
 Mₑ = ConstructM(lmax,s,m,c)
 F = eigen(Mₑ)
 vals = Complex.(F.values)::Array{Complex{Float64},1}
 vecs = F.vectors::Array{Complex{Float64},2}
 error, j = findmin(map(Anew -> abs2(Anew - Aₗₘprev), vals))
 vals[j], vecs[:,j], error
end
