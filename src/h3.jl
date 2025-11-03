#-----------------------------------------------------------------------------# H3Grid
struct H3Grid end

Base.show(io::IO, ::H3Grid) = print(io, styled"{bright_cyan: ⣿⣿ H3Grid}")

Base.getindex(g::H3Grid, i::Integer) = H3Cell(i, Int[])

function pentagons(g::H3Grid, res::Integer)
    out = Vector{UInt64}(undef, 12)
    LibH3.getPentagons(Cint(res), pointer(out))
    return H3Cell.(out)
end

#-----------------------------------------------------------------------------# H3Cell
"""
    H3Cell(index::UInt64)
    H3Cell(str::AbstractString)
    H3Cell(coord::AbstractCoordinates, res::Integer = 10)
"""
struct H3Cell
    index::UInt64
    function H3Cell(index::UInt64)
        LibH3.isValidCell(index) != 1 && throw(ArgumentError("Not a valid H3Cell Index: $index."))
        return new(index)
    end
end
function Base.show(io::IO, o::H3Cell)
    shape = is_pentagon(o) ? styled"{bright_red: ⬠}" : styled"{bright_green: ⬡}"
    print(io, styled"$shape{bright_cyan: H3Cell} {bright_magenta:res=$(resolution(o))} {bright_yellow:base_cell=$(base_cell(o))}")
    print(io, styled"{bright_black: $(digits(o))}")
end

Base.getindex(o::H3Cell, i::Integer) = H3Cell(base_cell(o), vcat(digits(o), i))

# constants from the H3 spec
const H3_NUM_DIGITS        = 15
const H3_DIGIT_BITS        = 3
const H3_DIGIT_MASK        = UInt64(0x7)     # 3 bits
const H3_BASE_CELL_BITS    = 7
const H3_RES_BITS          = 4
const H3_RESERVED3_BITS    = 3
const H3_MODE_BITS         = 4

# bit offsets from LSB (bit 0)
const H3_DIGITS_OFFSET     = 0
const H3_BASE_CELL_OFFSET  = H3_DIGITS_OFFSET + H3_NUM_DIGITS*H3_DIGIT_BITS      # 45
const H3_RES_OFFSET        = H3_BASE_CELL_OFFSET + H3_BASE_CELL_BITS             # 52
const H3_MODE_OFFSET       = H3_RES_OFFSET + H3_RES_BITS + H3_RESERVED3_BITS     # 59
const H3_MODE_CELL         = UInt64(1)

function H3Cell(base_cell::Integer, digits::Vector{<:Integer})
    0 ≤ base_cell ≤ 121 || throw(ArgumentError("base_cell must be in 0:121"))
    length(digits) ≤ H3_NUM_DIGITS || throw(ArgumentError("digits length must be ≤ $H3_NUM_DIGITS"))
    all(0 .≤ digits .≤ 6) || throw(ArgumentError("each digit must be in 0:6"))
    r = length(digits)
    bc = Int(base_cell)

    # pack digits (Digit 1 in the highest 3 bits of the 45-bit digit field)
    h = UInt64(0)
    @inbounds for i in 1:H3_NUM_DIGITS
        d = (i ≤ r) ? UInt64(digits[i]) : H3_DIGIT_MASK  # pad with 7
        shift = H3_DIGIT_BITS * (H3_NUM_DIGITS - i) + H3_DIGITS_OFFSET
        h |= (d & H3_DIGIT_MASK) << shift
    end

    # pack base cell, resolution, (reserved=0), mode
    h |= (UInt64(bc) & UInt64(0x7f)) << H3_BASE_CELL_OFFSET
    h |= (UInt64(r)  & UInt64(0x0f)) << H3_RES_OFFSET
    # reserved3 bits at RES_OFFSET + RES_BITS are left as 0
    h |= (H3_MODE_CELL  & UInt64(0x0f)) << H3_MODE_OFFSET

    return H3Cell(h)
end

function H3Cell(x::AbstractString)
    out = Ref{UInt64}()
    LibH3.stringToH3(String(x), out)
    return H3Cell(out[])
end

function H3Cell(coord::AbstractCoordinates, res::Integer = 10)
    out = Ref{UInt64}()
    LibH3.latLngToCell(Ref(LibH3.LatLng(coord)), Cint(res), out)
    return H3Cell(out[])
end

GI.geomtrait(::H3Cell) = GI.PolygonTrait()

function GI.centroid(::GI.PolygonTrait, o::H3Cell)
    g = Ref{LibH3.LatLng}()
    LibH3.cellToLatLng(o.index, g)
    return LonLat(LonLatAuthalic(g[]))
end

function GI.coordinates(::GI.PolygonTrait, o::H3Cell)
    g = Ref{LibH3.CellBoundary}()
    LibH3.cellToBoundary(o.index, g)
    out = [LonLat(LonLatAuthalic(x)) for x in g[].verts[1:g[].numVerts]]
    return vcat(out, out[1])  # close the polygon
end

GI.isgeometry(::H3Cell) = true
GI.crs(::H3Cell) = GFT.EPSG(4326)
GI.ncoord(::GI.PolygonTrait, o::H3Cell) = 2
GI.nhole(::GI.PolygonTrait, o::H3Cell) = 0
GI.ngeom(::GI.PolygonTrait, o::H3Cell) = 1
GI.getgeom(::GI.PolygonTrait, o::H3Cell, i::Integer) = GI.LineString(GI.coordinates(o))
GI.area(::GI.PolygonTrait, o::H3Cell) = (out = Ref{Cdouble}(); LibH3.cellAreaM2(o.index, out); out[])

#-----------------------------------------------------------------------------# H3Cell inspection
is_pentagon(o::H3Cell) = LibH3.isPentagon(o.index) == 1

# resolution(o::H3Cell) = Int(LibH3.getResolution(o.index))
@inline resolution(o::H3Cell) = Int((o.index >> H3_RES_OFFSET) & 0xF)  # Faster than using ccall

base_cell(o::H3Cell) = Int(LibH3.getBaseCellNumber(o.index))

h3string(o::H3Cell) = string(o.index, base=16)

@inline mode(o::H3Cell) = Int((o.index >> H3_MODE_OFFSET) & 0xF)

# Digit i is 1-based (Digit 1 is the highest 3 bits of the 45-bit digit field).
@inline function digit(h::UInt64, i::Int)
    @assert 1 ≤ i ≤ H3_NUM_DIGITS
    shift = H3_DIGIT_BITS * (H3_NUM_DIGITS - i)
    Int((h >> shift) & H3_DIGIT_MASK)  # 0–6 when used; 7 means "unused"
end

Base.digits(o::H3Cell) = [digit(o.index, i) for i in 1:resolution(o)]

# Faces that the cell intersects
function face_numbers(o::H3Cell)
    max_faces = Ref{Cint}()
    LibH3.maxFaceCount(o.index, max_faces)
    faces = Vector{Cint}(undef, max_faces[])
    LibH3.getIcosahedronFaces(o.index, pointer(faces))
    return filter!(!=(-1), faces)
end

function parent(o::H3Cell)
    digs = digits(o)
    isempty(digs) ? nothing :  H3Cell(base_cell(o), digs[1:end-1])
end

function children(o::H3Cell)
    size = Ref{Int64}()
    LibH3.cellToChildrenSize(o.index, Cint(resolution(o)+1), size)
    out = Vector{UInt64}(undef, size[])
    LibH3.cellToChildren(o.index, Cint(resolution(o)+1), pointer(out))
    return H3Cell.(out)
end

function siblings(o::H3Cell)
    digs = digits(o)
    isempty(digs) ? nothing : filter!(!=(o), children(parent(o)))
end

#-----------------------------------------------------------------------------# traversal
function grid_distance(a::H3Cell, b::H3Cell)
    dist = Ref{Int64}()
    LibH3.gridDistance(a.index, b.index, dist)
    return dist[]
end

haversine(a::H3Cell, b::H3Cell) = haversine(GI.centroid(a), GI.centroid(b))

destination(a::H3Cell, azimuth°, m) = H3Cell(destination(GI.centroid(a), azimuth°, m), resolution(a))

function grid_ring(o::H3Cell, k::Integer)
    k ≥ 0 || throw(ArgumentError("k must be non-negative.  Found: $k."))
    num_cells = Cint(6 * k)
    out = Vector{UInt64}(undef, num_cells)
    LibH3.gridRing(o.index, Cint(k), pointer(out))
    return [H3Cell(h) for h in out if h != 0]
end

function grid_disk(o::H3Cell, k::Integer)
    k ≥ 0 || throw(ArgumentError("k must be non-negative.  Found: $k."))
    num_cells = Cint(3 * k * (k + 1) + 1)
    out = Vector{UInt64}(undef, num_cells)
    LibH3.gridDisk(o.index, Cint(k), pointer(out))
    return [H3Cell(h) for h in out if h != 0]
end

function grid_path(a::H3Cell, b::H3Cell)
    size = Ref{Int64}()
    LibH3.gridPathCellsSize(a.index, b.index, size)
    out = Vector{UInt64}(undef, size[])
    LibH3.gridPathCells(a.index, b.index, pointer(out))
    filter!(!iszero, out)
    return H3Cell.(out)
end

#-----------------------------------------------------------------------------# cells
h3cells(geom, res::Integer = 10; kw...) = h3cells(GI.geomtrait(geom), geom, res; kw...)

h3cells(::GI.PointTrait, geom, res::Integer) = [H3Cell(LonLat(GI.coordinates(geom)), res)]

function h3cells(::GI.MultiPointTrait, geom, res::Integer)
    coords = LonLat.(GI.coordinates(geom))
    unique!(H3Cell.(coords, res))
end

function h3cells(::GI.LineTrait, geom, res::Integer; shortest_path = true)
    coords = LonLat.(GI.coordinates(geom))
    a = H3Cell(coords[1], res)
    b = H3Cell(coords[2], res)
    out = grid_path(a, b)
    if shortest_path
        return out
    else
        for cell in out, candidate in grid_ring(cell, 1)
            candidate in out && continue
            GO.disjoint(geom, candidate) && continue
            push!(out, candidate)
        end
        filter!(x -> !GO.disjoint(geom, x), out)
    end
    # # Now expand to include neighboring cells that the line crosses
    # if !shortest_path
    #     for cell in out, candidate in cell
    #         candidate in out && continue
    #         GO.disjoint(geom, candidate) && continue
    #         push!(out, candidate)
    #     end
    #     unique!(out)
    # end
    # out
end

function h3cells(::GI.LineStringTrait, geom, res::Integer; shortest_path=true)
    coords = GI.coordinates(geom)
    out = H3Cell[]
    @views for (a,b) in zip(coords[1:end-1], coords[2:end])
        line = GI.Line([a, b])
        union!(out, h3cells(line, res; shortest_path))
    end
    return out
end


# Get LibH3.GeoPolygon from GeoInterface polygon
function _h3polygon(geom)
    GI.trait(geom) == GI.PolygonTrait() || throw(ArgumentError("Expected Polygon geometry"))
    verts = map(GI.coordinates(geom)) do ring
        map(ring) do coord
            LibH3.LatLng(deg2rad(coord[2]), deg2rad(coord[1]))
        end
    end
    GC.@preserve verts begin
        geo_loops = LibH3.GeoLoop.(length.(verts), pointer.(verts))
    end
    GC.@preserve geo_loops begin
        n = length(geo_loops)
        return LibH3.GeoPolygon(geo_loops[1], n - 1, n > 1 ? pointer(geo_loops, 2) : C_NULL)
    end
end

function h3cells(::GI.PolygonTrait, geom, res::Integer; containment = :overlap)
    containment in (:center, :full, :overlap, :overlap_bbox) ||
        throw(ArgumentError("Invalid containment mode.  Expected one of `(:center, :full, :overlap, :overlap_bbox)`.  Found: $containment."))
    geo_polygon = _h3polygon(geom)
    flag_dict = Dict{Symbol, UInt32}(:center => 0, :full => 1, :overlap => 2, :overlap_bbox => 3)
    flag = flag_dict[containment]
    GC.@preserve geo_polygon begin
        max_n = Ref{Int64}()
        LibH3.maxPolygonToCellsSizeExperimental(Ref(geo_polygon), res, flag, max_n)
        out = zeros(UInt64, max_n[])
        LibH3.polygonToCellsExperimental(Ref(geo_polygon), res, flag, max_n[], out)
        H3Cell.(filter!(!iszero, unique!(out)))
    end
end

function h3cells(::GI.MultiPolygonTrait, geom, res::Integer; kw...)
    reduce(union, h3cells.(GI.getpolygon(geom), res; kw...))
end

function h3cells((; X, Y)::Extents.Extent, res::Integer; kw...)
    ls = GI.LineString([(X[1], Y[1]), (X[1], Y[2]), (X[2], Y[2]), (X[2], Y[1]), (X[1], Y[1])])
    h3cells(GI.Polygon([ls]), res; kw...)
end
