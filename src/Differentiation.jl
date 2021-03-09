struct LinearCombinationOf{T}
    values::Array{T,1}
    coeffs::Array{Complex{Float64},1}
end

function Base.(+)(x::HeunConfluentRadial,y::HeunConfluentRadial)
    LinearCombinationOf([x,y],[1,1])
end

function Base.(+)(x::HeunConfluentRadial,y::LinearCombinationOf{HeunConfluentRadial})
    if x ∈ y.values
        print("exists")
    else
        print("doesnt exist")
    end
end
