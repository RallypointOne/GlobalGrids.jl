using GlobalGrids
using Test

using GlobalGrids: IGEO7

@testset "GlobalGrids.jl" begin

@testset "IGEO7" begin
    x = IGEO7.encode(3, [1,2,3,4,5])
    (;base, digits) = IGEO7.decode(x)
    @test base == 2
    @test digits == [1:5..., fill(7,15)...]
end # IGEO7
end
