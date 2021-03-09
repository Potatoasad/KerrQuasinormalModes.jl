

"""
HeunConfluentLocalExpansion
    To Compute such a function one needs to decide the indices
    at each particlar point. There are 8 possible choices for this
    we must create an interface that allows people to choose which
    expansion they want
"""
struct HeunConfluentLocalExpansion
    p::Complex{Float64}; α::Complex{Float64}; γ::Complex{Float64};
    δ::Complex{Float64}; σ::Complex{Float64};
    acoeffs::Array{Complex{Float64},1}
    rcoeffs::Array{Complex{Float64},1}
end

function ComputeSeries(p,α,γ,δ,σ)
    αₖ(k) = -(k + 1)*(k + γ);
    βₖ(k) = k*(k - 4*p + γ + δ - 1);
    γₖ(k) = 4*p*(k + α - 1);
    an(k) = βₙ(k)/αₙ(k) ;
    bn(k) = γₙ(k)/αₙ(k);
    ComputeSeriesFromab(an,bn; N=250)
end
