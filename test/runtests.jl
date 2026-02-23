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

        pents = pentagons(H3Grid(), 0)
        for p in pents
            @test is_pentagon(p)
            @test length(GG.h3_face_numbers(p.index)) == 5
        end

        @test resolution(o) == res
        @test !is_pentagon(o)
        @test length(GG.digits(o)) == res
        for i in (res + 1):15
            @test GG.h3_digit(o.index, i) == 7
        end
        if res > 0
            p = GG.parent(o)
            @test resolution(p) == res - 1
            @test o in children(p)
            for sib in siblings(o)
                @test resolution(sib) == res
                @test GG.parent(sib) == p
            end
        else
            @test isnothing(GG.parent(o))
        end

        # Constructors
        @test o == H3Cell(encode(o))
        @test o == H3Cell(GI.centroid(o), resolution(o))
        @test o == H3Cell(o.index)
        @test o == H3Cell((coord[1], coord[2]), res)
        @test o == H3Cell(coord, res)

        # GeoInterface
        @test GI.area(o) > 0
        @test length(GI.coordinates(o)) == (is_pentagon(o) ? 6 : 7)
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
    @test haversine(c1, c2) > 0
    c3 = destination(c1, 90.0, 50_000.0)
    @test c3 isa H3Cell

    # base_cell and digits
    @test base_cell(c1) == GG.h3_base_cell(c1.index)
    @test GG.digits(c1) == GG.h3_digits(c1.index)
    @test encode(c1) == GG.h3_string(c1.index)

    # ncells
    @test ncells(H3Grid(), 0) == GG.h3_n_cells(0)
    @test ncells(H3Grid(), 5) == GG.h3_n_cells(5)

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

#-----------------------------------------------------------------------------# IGEO7Grid
@testset "IGEO7Grid" begin
    g = IGEO7Grid()
    @test g isa IGEO7Grid
    @test GI.isgeometry(g)
    @test GI.crs(g) == GI.crs(g)
    @test contains(sprint(show, g), "IGEO7Grid")
end

#-----------------------------------------------------------------------------# IGEO7Cell
@testset "IGEO7Cell" begin
    coord = LonLat(-75.0, 54.0)

    # Invalid resolution
    @test_throws ArgumentError IGEO7Cell(coord, -1)
    @test_throws ArgumentError IGEO7Cell(coord, 21)

    # Invalid base
    @test_throws ArgumentError IGEO7Cell(-1, Int[])
    @test_throws ArgumentError IGEO7Cell(12, Int[])

    # Invalid digit
    @test_throws ArgumentError IGEO7Cell(0, [7])

    # Pentagon cannot have digit 6
    @test_throws ArgumentError IGEO7Cell(0, [6])

    # Hex string constructor
    c = IGEO7Cell(coord, 3)
    @test IGEO7Cell(encode(c)) == c

    # Resolution sweep
    for res in 0:8
        o = IGEO7Cell(coord, res)

        @test resolution(o) == res
        # At res 0, all cells are pentagons (base vertices)
        if res == 0
            @test is_pentagon(o)
        else
            @test !is_pentagon(o)
        end
        @test length(GG.digits(o)) == res

        # Padding digits are 7
        for i in (res + 1):GG.IGEO7_NUM_DIGITS
            @test GG.igeo7_digit(o.index, i) == 7
        end

        # Parent/children
        if res > 0
            p = GG.parent(o)
            @test resolution(p) == res - 1
            @test o in children(p)
            for sib in siblings(o)
                @test resolution(sib) == res
                @test GG.parent(sib) == p
            end
        else
            @test isnothing(GG.parent(o))
        end

        # Round-trip: cell → centroid → cell
        @test o == IGEO7Cell(GI.centroid(o), res)

        # GeoInterface
        @test GI.area(o) > 0
        @test length(GI.coordinates(o)) == (is_pentagon(o) ? 6 : 7)
    end

    # Pentagon tests across resolutions
    for res in 0:5
        pents = pentagons(IGEO7Grid(), res)
        @test length(pents) == 12
        @test all(is_pentagon, pents)

        # Pentagon children count
        for p in pents
            if res < GG.IGEO7_MAX_RES
                ch = children(p)
                @test length(ch) == 6
                @test is_pentagon(ch[1])  # digit-0 child is pentagon
                @test all(!is_pentagon, ch[2:end])  # others are hexagons
            end
        end
    end

    # Hexagon children count
    hex = IGEO7Cell(0, [1])
    @test length(children(hex)) == 7

    # Cell count formula: 10*7^r + 2
    @test ncells(IGEO7Grid(), 0) == 12
    @test ncells(IGEO7Grid(), 1) == 72
    @test ncells(IGEO7Grid(), 2) == 492
    @test ncells(IGEO7Grid(), 3) == 3432

    # Area decreases with resolution
    areas = [GI.area(IGEO7Cell(coord, r)) for r in 0:5]
    @test issorted(areas, rev=true)
    @test all(>(0), areas)

    # IGEO7Cell haversine and destination
    c1 = IGEO7Cell(LonLat(0.0, 0.0), 5)
    c2 = IGEO7Cell(LonLat(1.0, 0.0), 5)
    @test haversine(c1, c2) > 0
    c3 = destination(c1, 90.0, 50_000.0)
    @test c3 isa IGEO7Cell

    # base_cell, digits, encode
    @test base_cell(c1) == GG.igeo7_base_cell(c1.index)
    @test GG.digits(c1) == GG.igeo7_digits(c1.index)
    @test encode(c1) == GG.igeo7_string(c1.index)

    # Show
    @test contains(sprint(show, c1), "IGEO7Cell")

    # decode
    @test decode(IGEO7Cell(0, [1, 2, 3])) == "0-123"

    # crosses_meridian
    @test !GG.crosses_meridian(IGEO7Cell(LonLat(0.0, 0.0), 5))
end

#-----------------------------------------------------------------------------# igeo7cells
@testset "igeo7cells" begin
    # Point
    x = igeo7cells(LonLat(0.0, 0.0), 5)
    @test length(x) == 1
    @test x[1] isa IGEO7Cell

    # MultiPoint
    x = igeo7cells(GI.MultiPoint([(0.0, 0.0), (10.0, 10.0), (20.0, 20.0)]), 3)
    @test length(x) >= 1
    @test all(c -> c isa IGEO7Cell, x)

    # Line
    x = igeo7cells(GI.Line([(0.0, 0.0), (5.0, 5.0)]), 3; containment=:center)
    @test length(x) >= 2
    @test all(c -> c isa IGEO7Cell, x)

    # LineString
    x = igeo7cells(GI.LineString([(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)]), 3)
    @test length(x) >= 2

    # Polygon (use large enough polygon for cells to have centroids inside)
    boundary = [(-20.0, -20.0), (-20.0, 20.0), (20.0, 20.0), (20.0, -20.0), (-20.0, -20.0)]
    poly = GI.Polygon([boundary])
    x = igeo7cells(poly, 2; containment=:center)
    @test length(x) >= 1
    @test all(c -> c isa IGEO7Cell, x)

    # Invalid containment
    @test_throws ArgumentError igeo7cells(GI.Line([(0.0, 0.0), (1.0, 1.0)]), 3; containment=:invalid)
    @test_throws ArgumentError igeo7cells(poly, 3; containment=:invalid)
end

#-----------------------------------------------------------------------------# DGGGrid
@testset "DGGGrid" begin
    g = ISEA3H()
    @test g isa ISEA3H
    @test g isa DGGGrid{:isea, 3, :hex}
    @test GI.isgeometry(g)
    @test GI.crs(g) == GI.crs(g)
    @test contains(sprint(show, g), "DGGGrid")
end

#-----------------------------------------------------------------------------# DGGCell
@testset "DGGCell" begin
    coord = LonLat(-75.0, 54.0)

    # Invalid resolution
    @test_throws ArgumentError ISEA3HCell(coord, -1)
    @test_throws ArgumentError ISEA3HCell(coord, 29)

    # Invalid base
    @test_throws ArgumentError ISEA3HCell(-1, Int[])
    @test_throws ArgumentError ISEA3HCell(12, Int[])

    # Invalid digit
    @test_throws ArgumentError ISEA3HCell(0, [3])  # aperture-3: digits must be 0–2

    # Hex string constructor
    c = ISEA3HCell(coord, 3)
    @test ISEA3HCell(encode(c)) == c

    # Tuple constructor
    @test ISEA3HCell((coord[1], coord[2]), 3) == c

    # Resolution sweep
    for res in 0:12
        o = ISEA3HCell(coord, res)

        @test resolution(o) == res
        if res == 0
            @test is_pentagon(o)
        else
            @test !is_pentagon(o)
        end
        @test length(GG.digits(o)) == res

        # Padding digits are 0x3 (all-1s mask for 2-bit digits)
        for i in (res + 1):GG._dgg_maxd(Val(3))
            @test GG.dgg_digit(o.index, i, Val(3)) == 3
        end

        # Parent/children
        if res > 0
            p = GG.parent(o)
            @test resolution(p) == res - 1
            @test o in children(p)
            for sib in siblings(o)
                @test resolution(sib) == res
                @test GG.parent(sib) == p
            end
        else
            @test isnothing(GG.parent(o))
        end

        # Round-trip: cell → centroid → cell
        @test o == ISEA3HCell(GI.centroid(o), res)

        # GeoInterface
        @test GI.area(o) > 0
        @test length(GI.coordinates(o)) == (is_pentagon(o) ? 6 : 7)
    end

    # 12 pentagons at each resolution
    for res in 0:5
        pents = pentagons(ISEA3H(), res)
        @test length(pents) == 12
        @test all(is_pentagon, pents)
    end

    # Cell count formula: 10*3^r + 2
    @test ncells(ISEA3H(), 0) == 12
    @test ncells(ISEA3H(), 1) == 32
    @test ncells(ISEA3H(), 2) == 92
    @test ncells(ISEA3H(), 3) == 272

    # Area decreases with resolution
    areas = [GI.area(ISEA3HCell(coord, r)) for r in 0:5]
    @test issorted(areas, rev=true)
    @test all(>(0), areas)

    # DGGCell haversine and destination (use res 8 so cells are small enough)
    c1 = ISEA3HCell(LonLat(0.0, 0.0), 8)
    c2 = ISEA3HCell(LonLat(1.0, 0.0), 8)
    @test haversine(c1, c2) > 0
    c3 = destination(c1, 90.0, 50_000.0)
    @test c3 isa ISEA3HCell

    # base_cell, digits, encode
    @test base_cell(c1) == GG.dgg_base(c1.index)
    @test GG.digits(c1) == GG.dgg_digits(c1.index, Val(3))
    @test encode(c1) == GG.dgg_string(c1.index)

    # Show
    @test contains(sprint(show, c1), "DGGCell")

    # Decode
    @test decode(ISEA3HCell(0, [1, 2, 0])) == "0-120"

    # Round-trip with a grid of test points at multiple resolutions
    test_lons = range(-180, 180, length=20)
    test_lats = range(-85, 85, length=10)
    for res in [1, 3, 5, 8]
        for lon in test_lons, lat in test_lats
            c = ISEA3HCell(LonLat(Float64(lon), Float64(lat)), res)
            @test c == ISEA3HCell(GI.centroid(c), res)
        end
    end
end

#-----------------------------------------------------------------------------# dggcells
@testset "dggcells" begin
    # Point
    x = dggcells(LonLat(0.0, 0.0), 5)
    @test length(x) == 1
    @test x[1] isa ISEA3HCell

    # MultiPoint
    x = dggcells(GI.MultiPoint([(0.0, 0.0), (10.0, 10.0), (20.0, 20.0)]), 3)
    @test length(x) >= 1
    @test all(c -> c isa ISEA3HCell, x)

    # Line
    x = dggcells(GI.Line([(0.0, 0.0), (5.0, 5.0)]), 3; containment=:center)
    @test length(x) >= 2
    @test all(c -> c isa ISEA3HCell, x)

    # LineString
    x = dggcells(GI.LineString([(0.0, 0.0), (5.0, 5.0), (10.0, 0.0)]), 3)
    @test length(x) >= 2

    # Polygon (use large enough polygon for cells to have centroids inside)
    boundary = [(-20.0, -20.0), (-20.0, 20.0), (20.0, 20.0), (20.0, -20.0), (-20.0, -20.0)]
    poly = GI.Polygon([boundary])
    x = dggcells(poly, 2; containment=:center)
    @test length(x) >= 1
    @test all(c -> c isa ISEA3HCell, x)

    # Invalid containment
    @test_throws ArgumentError dggcells(GI.Line([(0.0, 0.0), (1.0, 1.0)]), 3; containment=:invalid)
    @test_throws ArgumentError dggcells(poly, 3; containment=:invalid)
end
