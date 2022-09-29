module KerrQuasinormalModes

using LinearAlgebra
using NLsolve
using Statistics
using CSV
using StaticArrays

using Parameters
using BSplineKit

#=
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
=#


include(joinpath(@__DIR__, "ModeSolver" ,"Angular.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"Radial.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"RootSolver.jl"))
include(joinpath(@__DIR__, "ModeSolver/Schwarszchild" ,"Schwarszchild.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"SpinSequenceOptions.jl"))
include(joinpath(@__DIR__, "ModeSolver" ,"SpinSequence.jl"))


include(joinpath(@__DIR__, "ModeFunctionInterface" ,"SpinWeightedSphericalLookup.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"Interface.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"ConfluentHeun.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"LinearCombinations.jl"))
include(joinpath(@__DIR__, "ModeFunctionInterface" ,"Derivative.jl"))

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
export CP
export ModeSequence

#For debugging
export HeunConfluentLocalExpansion, ComputeSeries, qnm, Custom
#export importqnm, qnmfunctionnew
end
