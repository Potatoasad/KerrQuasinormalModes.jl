module KerrQuasinormalModes

using LinearAlgebra
using NLsolve
using Statistics
using CSV
using StaticArrays

include("AngularFunctions.jl")
include("ParameterFunctions.jl")
include("ContinuedFractionFunctions.jl")
include("SolverFunctions.jl")
include("Schwarzschild.jl")
include("ModeCalculator.jl")
include("QnmRotationSeries.jl")
include("SpinWeightedSphericalLookup.jl")
include("Interface.jl")
include("ConfluentHeun.jl")
include("LinearCombinations.jl")
include("Derivative.jl")

#Ψ = qnmfunction(s=-2,l=2,n=2,m=-2,a=0.4)
#Ψ(2.0)
function asad()
    print("blah blah blo new yesterday heres somthing new")
end

export SpinWeightedSpherical
export HeunConfluentRadial
export qnmfunction
export QuasinormalModeFunction
export asad
export ComplexPlot

#For debugging
export HeunConfluentLocalExpansion, ComputeSeries, qnm, Custom
export importqnm, qnmfunctionnew
end
