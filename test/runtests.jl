using KerrQuasinormalModes
using Test

function TestIfValueOfRadialWorks()
    a = -0.402546147858408 - 0.06278653515072422*im;
    Ψ = qnmfunction(s=-2,l=2,n=2,m=-2,a=0.4)
    (abs(Ψ(2.0)-a) < 10^(-9))
end

function TestIfQnmExecutes()
    Ψ = qnmfunction(s=-2,l=2,n=2,m=-2,a=0.4)
    r = 2.0:0.01:6.0
    Ψ.(r)
    true
end

@testset "QuasinormalModes.jl" begin
    # Write your tests here.
    @test TestIfValueOfRadialWorks()
    @test TestIfQnmExecutes()
end
