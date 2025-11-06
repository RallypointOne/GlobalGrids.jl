using GlobalGrids
using Test

import GlobalGrids as GG
import GeometryBasics as GB
import GeoInterface as GI
import GeometryOps as GO


#-----------------------------------------------------------------------------# H3Grid
@testset "H3Grid" begin
    g = H3Grid()
    for i in 0:121
        @test g[i] isa H3Cell
    end
    @test_throws Exception g[-1]
    @test_throws Exception g[122]
end

#-----------------------------------------------------------------------------# H3Cell constructors
@testset "H3Cell" begin
    coord = LonLat(-75.0, 54.0)
    @test_throws Exception H3Cell(coord, -1)
    @test_throws Exception H3Cell(coord, 16)

    for res in 0:15
        o = H3Cell(coord, res)

        pentagons = GG.h3_pentagons(0)
        for p in pentagons
            @test GG.is_pentagon(p)
            @test length(GG.h3_face_numbers(p.index)) == 5
        end

        @test GG.resolution(o) == res
        @test !GG.is_pentagon(o)
        @test length(GG.h3_digits(o.index)) == res
        for i in (res + 1):15
            @test GG.h3_digit(o.index, i) == 7
        end
        if res > 0
            p = GG.parent(o)
            @test GG.resolution(p) == res - 1
            @test o in GG.children(p)
            for sib in GG.siblings(o)
                @test GG.resolution(sib) == res
                @test GG.parent(sib) == p
            end
        else
            @test isnothing(GG.parent(o))
        end

        # Constructors
        @test o == H3Cell(GG.h3_string(o))
        @test o == H3Cell(GI.centroid(o), GG.resolution(o))
        @test o == H3Cell(o.index)
        @test o == H3Cell((coord[1], coord[2]), res)
        @test o == H3Cell(coord, res)

        # GeoInterface
        @test GI.area(o) > 0
        @test length(GI.coordinates(o)) == 7  # closed ring
        @test GO.contains(o, GI.centroid(o))

        # operations
        @test GG.grid_distance(o, o) == 0
        @test length(GG.grid_ring(o, 1)) == 6
        @test length(GG.grid_disk(o, 1)) == 7

        if res > 5
            for k in 0:5
                for cell in GG.grid_ring(o, k)
                    @test GG.grid_distance(o, cell) == k
                    @test length(GG.grid_path(o, cell)) == k + 1
                end
            end
        end
    end
end

#-----------------------------------------------------------------------------# h3cells
@testset "h3cells" begin
    boundary = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]
    inner = [(0.2, 0.2), (0.2, 0.8), (0.8, 0.8), (0.8, 0.2), (0.2, 0.2)]
    polygon = GI.Polygon([boundary, inner])
    multipolygon = GI.MultiPolygon([polygon, GI.Polygon([[(5.0, 5.0), (5.0, 6.0), (6.0, 6.0), (6.0, 5.0), (5.0, 5.0)]])])

    # Point
    x = h3cells(boundary[1], 5)
    @test GO.contains(x[1], boundary[1])

    # Multipoint
    x = h3cells(GI.MultiPoint(boundary[1:3]), 5)

    # Line
    x = h3cells(GI.Line(boundary[1:2]), 5; containment = :shortest_path)
    x = h3cells(GI.Line(boundary[1:2]), 5; containment = :overlap)

    # Linestring
    x = h3cells(GI.LineString(boundary), 5; containment = :shortest_path)

    # Polygon
    x = h3cells(polygon, 5; containment = :center)

    # Multipolygon
    x = h3cells(multipolygon, 5; containment = :overlap)

    # Extent
    x = h3cells(GI.extent(polygon), 2; containment = :overlap_bbox)

end

# @testset "GlobalGrids.jl" begin

# #-----------------------------------------------------------------------------# Coordinates
# @testset "Coordinates" begin
#     origin = GG.LonLat(-75.0, 54.0)
#     types = [GG.LonLat, GG.LonLatAuthalic, GG.ECEF, GG.ISEA, GG.ISEACube]
#     for T in types, S in types
#         # @info "Testing transform: $(T) -> $(S) -> $(T)"
#         x = T(origin)
#         y = S(x)
#         x2 = T(y)
#         @test x.pt ≈ x2.pt
#     end
# end # Coordinates

# #-----------------------------------------------------------------------------# H3
# @testset "H3" begin
#     # Sanity checks:
#     g = GG.H3Grid()
#     for res in 0:15
#         @test length(GG.pentagons(g, res)) == 12
#     end
#     @test g isa GG.H3Grid
#     ll = GG.LonLat(-75.0, 54.0)

#     # resolution 0
#     o = GG.H3Cell(ll, 0)
#     @test o isa GG.H3Cell
#     @test GG.resolution(o) == 0
#     @test isnothing(GG.parent(o))
#     @test isnothing(GG.siblings(o))

#     for res in 1:15
#         o = GG.H3Cell(ll, res)
#         o2 = GG.H3Cell(GG.h3string(o))
#         @test GI.area(o) > 0
#         @test o == o2
#         @test o isa GG.H3Cell
#         @test GG.resolution(o) == res
#         @test !GG.is_pentagon(o)
#         for sib in GG.siblings(o)
#             @test GG.grid_distance(o, sib) in (1, 2)
#             @test length(GG.grid_path(o, sib)) == GG.grid_distance(o, sib) + 1
#             @test GG.haversine(o, sib) > 0.0
#         end
#         for cell in GG.grid_ring(o, 1)
#             @test GG.grid_distance(o, cell)  == 1
#         end
#         for cell in GG.grid_disk(o, 2)
#             @test GG.grid_distance(o, cell) ≤ 2
#         end
#     end
#     @testset "h3cells" begin
#         # Point
#         ll = GG.LonLat(-75.0, 54.0)
#         for res in 0:15
#             @test GG.h3cells(ll, res) isa Vector{GG.H3Cell}
#         end

#         # MultiPoint/Line
#         mll = GB.MultiPoint{2, Float64}([GB.Point2(-75.0, 54.0), GB.Point2(-80.0, 50.0)])
#         line = GI.Line(GI.coordinates(mll))
#         for res in 0:15
#             @test GG.h3cells(mll, res) isa Vector{GG.H3Cell}
#             @test GG.h3cells(line, res) isa Vector{GG.H3Cell}
#         end

#         # LineString
#         linestring = GI.LineString([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (0.0, 0.0)])
#         for res in 0:15
#             @test GG.h3cells(linestring, res) isa Vector{GG.H3Cell}
#         end

#         # Polygon
#         o = GG.H3Cell(ll)
#         ring = GG.grid_ring(o, 3)
#         vertices = [GI.centroid(c) for c in ring]
#         poly = GI.Polygon([[vertices..., vertices[1]]])
#         for res in 0:15
#             @test length(GG.h3cells(poly, res)) > 0
#         end

#         # MultiPolygon
#         vertices2 = [x .+ GG.LonLat(5.0, 5.0) for x in vertices]
#         poly2 = GI.Polygon([[vertices2..., vertices2[1]]])
#         multipoly = GI.MultiPolygon([poly, poly2])
#         for res in 0:15
#             @test length(GG.h3cells(multipoly, res)) > 0
#         end

#         # Extents
#         ex = GI.extent(poly)
#         for res in 0:15
#             @test length(GG.h3cells(ex, res)) > 0
#         end

#     end  # cells
# end  # H3


# end # GlobalGrids.jl
