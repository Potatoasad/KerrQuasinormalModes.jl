using CSV

SpinLookUp = Dict{Tuple{Int64,Int64,Int64},String}()

for row in CSV.File(joinpath(@__DIR__, "SpinWeightedSphericalHarmonicsTable.csv"))
    s = row.s; l=row.l; m=row.m; func = row.Func
    SpinLookUp[(s,l,m)] = row.Func
end

function sterlings(n)
    return sqrt(2*pi*n)*(n/exp(1))^n
end

function MakeSwitchCase(param,iterable,args)
    newstr = """if """
    for x ∈ iterable
       newstr *= """$(param) == $(x)
    $(getcode(args...,x))
elseif """
    end
    newstr[1:(end-7)]*" end"
end


function getcode(s,l,m)
    return "return begin $(SpinLookUp[(s,l,m)]) end"
end

getcode(s,l) = MakeSwitchCase("m",-l:l,(s,l))
getcode(s) = MakeSwitchCase("l",abs(s):20,(s,))
getcode() = MakeSwitchCase("s",-2:2,())

#print("again")
"""
function SpinWeightedSphericalCalculation(z,s,l,m)
    zm = 1-z; zp = 1+z;
    if l < max(abs(s),abs(m))
        return 0.0
    end
    #println("BEFORE WE HAVE",s," ", l," ", m)
    $(getcode())
    if ((l+m > 20) || (l+s > 20) || (l-m > 20) || (l-s > 20))
        term1 = (2*l+1)*(sterlings(l+m)/sterlings(l+s))*(sterlings(l-m)/sterlings(l-s))
    else
        term1 = (2*l+1)*(factorial(l+m)/factorial(l+s))*(factorial(l-m)/factorial(l-s))
    end
    println("HERE WE HAVE",s," ", l," ", m)
    A = sqrt(term1/(4*pi))*(-1)^m

    #The sin term right outside the summation
    sinterm = ((1-z)/2)^l
    #The summation terms
    sumterms = 0.0
    for r = 0:(l-s)
        consts = binomial(l-s,r)*binomial(l+s,r+s-m)*(-1)^(l-r-s)
        cotterm = (sqrt((1+z)/(1-z)))^((2*r+s-m))
        sumterms += consts*cotterm
    end
    return A*sinterm*sumterms
end
""" |> Meta.parse |> eval
