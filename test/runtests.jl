using GlobalGrids
using Test

import GlobalGrids as GG
import GeometryBasics as GB
import GeoInterface as GI
import GeometryOps as GO


#-----------------------------------------------------------------------------# icosahedron
@testset "icosahedron" begin
    ico = icosahedron()
    @test ico isa GB.Mesh
    @test length(ico.position) == 12
    @test length(GB.faces(ico)) == 20

    # Preset orientations
    for preset in [:poles, :isea, :dymaxion]
        m = icosahedron(preset)
        @test m isa GB.Mesh
        @test length(m.position) == 12
    end
    @test_throws ArgumentError icosahedron(:invalid)

    # Custom vertex/azimuth
    ico2 = icosahedron(GB.Point3d(0, 0, -1), 0.0)
    @test ico2 isa GB.Mesh

    # Dual
    d = GG.dual(ico)
    @test d isa GB.Mesh
    @test length(d.position) == 20
end

#-----------------------------------------------------------------------------# LonLat
@testset "LonLat" begin
    p = LonLat(-75.0, 54.0)
    @test p.lon == -75.0
    @test p.lat == 54.0
    @test p[1] == -75.0
    @test p[2] == 54.0
    @test_throws BoundsError p[3]
    @test length(p) == 2
    @test eltype(p) == Float64
    @test collect(p) == [-75.0, 54.0]

    # Tuple constructor
    p2 = LonLat((1.0, 2.0))
    @test p2.lon == 1.0 && p2.lat == 2.0

    # GeoInterface
    @test GI.isgeometry(p)
    @test GI.geomtrait(p) == GI.PointTrait()
    @test GI.ncoord(p) == 2
    @test GI.getcoord(p, 1) == -75.0
    @test GI.getcoord(p, 2) == 54.0

    # Show
    @test contains(sprint(show, p), "LonLat")
end

#-----------------------------------------------------------------------------# haversine, destination, azimuth
@testset "geodesic functions" begin
    a = LonLat(0.0, 0.0)
    b = LonLat(0.0, 1.0)

    # haversine: ~111 km per degree of latitude
    d = GG.haversine(a, b)
    @test 110_000 < d < 112_000

    # Same point → zero distance
    @test GG.haversine(a, a) == 0.0

    # Azimuth: due north
    @test GG.azimuth(a, b) ≈ 0.0 atol=1e-10

    # Azimuth: due east
    c = LonLat(1.0, 0.0)
    @test GG.azimuth(a, c) ≈ 90.0 atol=0.1

    # Destination: go north 111 km ≈ 1 degree of latitude
    dest = GG.destination(a, 0.0, 111_195.0)
    @test dest.lat ≈ 1.0 atol=0.01
    @test dest.lon ≈ 0.0 atol=0.01
end

#-----------------------------------------------------------------------------# crosses_meridian
@testset "crosses_meridian" begin
    # Cell near prime meridian — should not cross 180
    c = H3Cell(LonLat(0.0, 0.0), 5)
    @test !GG.crosses_meridian(c)
end

#-----------------------------------------------------------------------------# H3Grid
@testset "H3Grid" begin
    g = H3Grid()
    @test g isa H3Grid
    @test GI.isgeometry(g)
    @test GI.crs(g) == GI.crs(g)

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
        @test length(GI.coordinates(o)) == (GG.is_pentagon(o) ? 6 : 7)
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

    # H3Cell haversine and destination
    c1 = H3Cell(LonLat(0.0, 0.0), 5)
    c2 = H3Cell(LonLat(1.0, 0.0), 5)
    @test GG.haversine(c1, c2) > 0
    c3 = GG.destination(c1, 90.0, 50_000.0)
    @test c3 isa H3Cell

    # Show
    @test contains(sprint(show, c1), "H3Cell")
end

#-----------------------------------------------------------------------------# h3cells
@testset "h3cells" begin
    boundary = [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]
    inner = [(0.2, 0.2), (0.2, 0.8), (0.8, 0.8), (0.8, 0.2), (0.2, 0.2)]
    polygon = GI.Polygon([boundary, inner])
    multipolygon = GI.MultiPolygon([polygon, GI.Polygon([[(5.0, 5.0), (5.0, 6.0), (6.0, 6.0), (6.0, 5.0), (5.0, 5.0)]])])

    # Point
    x = h3cells(boundary[1], 5)
    @test length(x) == 1
    @test x[1] isa H3Cell
    @test GO.contains(x[1], boundary[1])

    # Multipoint
    x = h3cells(GI.MultiPoint(boundary[1:3]), 5)
    @test length(x) >= 1
    @test all(c -> c isa H3Cell, x)

    # Line
    x = h3cells(GI.Line(boundary[1:2]), 5; containment = :shortest_path)
    @test length(x) >= 2
    x = h3cells(GI.Line(boundary[1:2]), 5; containment = :overlap)
    @test length(x) >= 1

    # Linestring
    x = h3cells(GI.LineString(boundary), 5; containment = :shortest_path)
    @test length(x) >= 4

    # Polygon
    x = h3cells(polygon, 5; containment = :center)
    @test length(x) >= 1
    @test all(c -> c isa H3Cell, x)

    # Multipolygon
    x = h3cells(multipolygon, 5; containment = :overlap)
    @test length(x) >= 1

    # Extent
    x = h3cells(GI.extent(polygon), 2; containment = :overlap_bbox)
    @test length(x) >= 1

    # Invalid containment
    @test_throws ArgumentError h3cells(GI.Line(boundary[1:2]), 5; containment = :invalid)
    @test_throws ArgumentError h3cells(polygon, 5; containment = :invalid)
end
