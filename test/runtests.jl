using QuasinormalModes
using Test

function TestIfValueOfRadialWorks()
    a = -0.402546147858408 - 0.06278653515072422*im;
    Ψ = qnmfunction(s=-2,l=2,n=2,m=-2,a=0.4)
    (abs(Ψ(2.0)-a) < 10^(-9))
end

@testset "QuasinormalModes.jl" begin
    # Write your tests here.
    @test TestIfValueOfRadialWorks()
end
