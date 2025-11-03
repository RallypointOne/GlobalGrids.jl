using GlobalGrids
using Test

import GlobalGrids as GG
import GeometryBasics as GB
import GeoInterface as GI

@testset "GlobalGrids.jl" begin

#-----------------------------------------------------------------------------# Coordinates
@testset "Coordinates" begin
    origin = GG.LonLat(-75.0, 54.0)
    types = [GG.LonLat, GG.LonLatAuthalic, GG.ECEF, GG.ISEA, GG.ISEACube]
    for T in types, S in types
        # @info "Testing transform: $(T) -> $(S) -> $(T)"
        x = T(origin)
        y = S(x)
        x2 = T(y)
        @test x.pt ≈ x2.pt
    end
end # Coordinates

#-----------------------------------------------------------------------------# H3
@testset "H3" begin
    # Sanity checks:
    g = GG.H3Grid()
    for res in 0:15
        @test length(GG.pentagons(g, res)) == 12
    end
    @test g isa GG.H3Grid
    ll = GG.LonLat(-75.0, 54.0)

    # resolution 0
    o = GG.H3Cell(ll, 0)
    @test o isa GG.H3Cell
    @test GG.resolution(o) == 0
    @test isnothing(GG.parent(o))
    @test isnothing(GG.siblings(o))

    for res in 1:15
        o = GG.H3Cell(ll, res)
        o2 = GG.H3Cell(GG.h3string(o))
        @test GI.area(o) > 0
        @test o == o2
        @test o isa GG.H3Cell
        @test GG.resolution(o) == res
        @test !GG.is_pentagon(o)
        for sib in GG.siblings(o)
            @test GG.grid_distance(o, sib) in (1, 2)
            @test length(GG.grid_path(o, sib)) == GG.grid_distance(o, sib) + 1
            @test GG.haversine(o, sib) > 0.0
        end
        for cell in GG.grid_ring(o, 1)
            @test GG.grid_distance(o, cell)  == 1
        end
        for cell in GG.grid_disk(o, 2)
            @test GG.grid_distance(o, cell) ≤ 2
        end
    end
    @testset "h3cells" begin
        # Point
        ll = GG.LonLat(-75.0, 54.0)
        for res in 0:15
            @test GG.h3cells(ll, res) isa Vector{GG.H3Cell}
        end

        # MultiPoint/Line
        mll = GB.MultiPoint{2, Float64}([GB.Point2(-75.0, 54.0), GB.Point2(-80.0, 50.0)])
        line = GI.Line(GI.coordinates(mll))
        for res in 0:15
            @test GG.h3cells(mll, res) isa Vector{GG.H3Cell}
            @test GG.h3cells(line, res) isa Vector{GG.H3Cell}
        end

        # LineString
        linestring = GI.LineString([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.0)])
        for res in 0:15
            @test GG.h3cells(linestring, res) isa Vector{GG.H3Cell}
        end

        # Polygon
        o = GG.H3Cell(ll)
        ring = GG.grid_ring(o, 3)
        vertices = [GI.centroid(c) for c in ring]
        poly = GI.Polygon([[vertices..., vertices[1]]])
        for res in 0:15
            @test length(GG.h3cells(poly, res)) > 0
        end

        # MultiPolygon
        vertices2 = [x .+ GG.LonLat(5.0, 5.0) for x in vertices]
        poly2 = GI.Polygon([[vertices2..., vertices2[1]]])
        multipoly = GI.MultiPolygon([poly, poly2])
        for res in 0:15
            @test length(GG.h3cells(multipoly, res)) > 0
        end

        # Extents
        ex = GI.extent(poly)
        for res in 0:15
            @test length(GG.h3cells(ex, res)) > 0
        end

    end  # cells
end  # H3


end # GlobalGrids.jl
