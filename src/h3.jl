#-----------------------------------------------------------------------------# utils for libh3
const H3Error = UInt32
const H3Index = UInt64

const H3_N_PENTAGONS        = 12
const H3_N_BASE_CELLS       = 122
const H3_NUM_DIGITS         = 15
const H3_DIGIT_BITS         = 3
const H3_DIGIT_MASK         = UInt64(0x7)  # 3 bits
const H3_BASE_CELL_BITS     = 7
const H3_RES_BITS           = 4
const H3_RESERVED3_BITS     = 3
const H3_MODE_BITS          = 4

# bit offsets from LSB (bit 0)
const H3_DIGITS_OFFSET      = 0
const H3_BASE_CELL_OFFSET   = H3_DIGITS_OFFSET + H3_NUM_DIGITS*H3_DIGIT_BITS    # 45
const H3_RES_OFFSET         = H3_BASE_CELL_OFFSET + H3_BASE_CELL_BITS           # 52
const H3_MODE_OFFSET        = H3_RES_OFFSET + H3_RES_BITS + H3_RESERVED3_BITS   # 59
const H3_MODE_CELL          = UInt64(1)

const H3_PENTAGON_BASE_CELLS = SVector((4, 14, 24, 38, 49, 58, 63, 72, 83, 97, 107, 117))

struct H3LatLng
    lat::Cdouble
    lng::Cdouble
end
LonLat(o::H3LatLng) = LonLat(rad2deg(o.lng), rad2deg(o.lat))
H3LatLng(o::LonLat) = H3LatLng(deg2rad(o.lat), deg2rad(o.lon))

struct H3CellBoundary
    numVerts::Cint
    verts::NTuple{10, H3LatLng}
end

# Note: verts[1] does not need to be repeated at verts[end]
struct H3GeoLoop
    numVerts::Cint
    verts::Ptr{H3LatLng}
end

struct H3GeoPolygon
    geoloop::H3GeoLoop
    numHoles::Cint
    holes::Ptr{H3GeoLoop}
end


function h3check(err::H3Error)
    if err != 0
        str = unsafe_string(@ccall(libh3.describeH3Error(err::H3Error)::Ptr{Cchar}))
        throw(ErrorException(styled"libh3 call failed with error code $err: {bright_red:$str}"))
    end
end

#-----------------------------------------------------------------------------# H3Grid
"""
    H3Grid()

The H3 hexagonal discrete global grid system developed by Uber.  Uses CRS
EPSG:4326 (WGS84 longitude/latitude).
"""
struct H3Grid <: AbstractGrid end

GI.crs(::H3Grid) = GFT.EPSG(4326)
GI.ncoord(::H3Grid) = 2
GI.npolygon(::GI.PolyhedralSurfaceTrait, o::H3Grid) = h3_n_cells(h3_resolution(o.index))
GI.ngeom(::GI.PolyhedralSurfaceTrait, grid::H3Grid) = GI.npolygon(grid)

#-----------------------------------------------------------------------------# H3Cell
"""
    H3Cell(index::UInt64)
    H3Cell(base::Integer, digits::AbstractVector{<:Integer})
    H3Cell(coord::LonLat, res::Integer=10)
    H3Cell(str::AbstractString)

A single cell in the H3 grid.  Wraps a `UInt64` index that encodes the base
cell (0–121), resolution (0–15), and refinement digits (0–6).

### Examples

    H3Cell(LonLat(-75.0, 54.0), 7)
    H3Cell(121, [1, 5, 0])
"""
struct H3Cell <: AbstractCell
    index::UInt64
    function H3Cell(idx::UInt64; validate::Bool = true)
        if validate
            @ccall(libh3.isValidCell(idx::UInt64)::Cint) == 1 || throw(ArgumentError("Not a valid H3Cell Index: $idx."))
        end
        new(idx)
    end
end

Base.getindex(g::H3Grid, i::Integer) = H3Cell(i, Int[])
Base.getindex(c::H3Cell, i::Integer) = H3Cell(h3_base_cell(c.index), [h3_digits(c.index)..., i])

function H3Cell(base::Integer, digits::AbstractVector{<:Integer})
    0 ≤ base ≤ 121 || throw(ArgumentError("base must be in 0:121"))
    length(digits) ≤ H3_NUM_DIGITS || throw(ArgumentError("digits length must be ≤ $H3_NUM_DIGITS"))
    all(x -> 0 ≤ x ≤ 6, digits) || throw(ArgumentError("each digit must be in 0:6"))
    r = length(digits)
    bc = Int(base)
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
    H3Cell(h)
end

function H3Cell(coord::H3LatLng, res::Integer = 10)
    out = Ref{UInt64}()
    h3check(@ccall libh3.latLngToCell(Ref(coord)::Ptr{H3LatLng}, Cint(res)::Cint, out::Ptr{UInt64})::H3Error)
    return H3Cell(out[])
end

H3Cell(coord::LonLat, res::Integer=10) = H3Cell(H3LatLng(coord), res)

H3Cell((x, y)::NTuple{2, Real}, res::Integer = 10) = H3Cell(LonLat(x, y), res)

H3Cell(str::AbstractString) = H3Cell(parse(UInt64, str; base=16))

#-----------------------------------------------------------------------------# GeoInterface
function GI.centroid(::GI.PolygonTrait, o::H3Cell)
    g = Ref{H3LatLng}()
    h3check(@ccall libh3.cellToLatLng(o.index::H3Index, g::Ptr{H3LatLng})::H3Error)
    return LonLat(g[])
end
function GI.coordinates(::GI.PolygonTrait, o::H3Cell)
    g = Ref{H3CellBoundary}()
    h3check(@ccall libh3.cellToBoundary(o.index::H3Index, g::Ptr{H3CellBoundary})::H3Error)
    out = [LonLat(x) for x in g[].verts[1:g[].numVerts]]
    return vcat(out, out[1])  # close the polygon
end
function GI.area(::GI.PolygonTrait, o::H3Cell)
    out = Ref{Cdouble}()
    h3check(@ccall(libh3.cellAreaM2(o.index::UInt64, out::Ptr{Cdouble})::H3Error))
    return out[]
end

#-----------------------------------------------------------------------------# AbstractCell interface
grid(o::H3Cell) = H3Grid()

# resolution
h3_resolution(x::UInt64) = Int((x >> H3_RES_OFFSET) & 0xF)  # Faster than using ccall
resolution(o::H3Cell) = h3_resolution(o.index)

# icon
h3_is_pentagon(x::UInt64) = @ccall(libh3.isPentagon(x::UInt64)::Cint) == 1
is_pentagon(o::H3Cell) = h3_is_pentagon(o.index)
icon(o::H3Cell) = is_pentagon(o) ? styled"{bright_red: ⬠}" : styled"{bright_green: ⬡}"

# decode
h3_base_cell(x::UInt64) = Int((x >> H3_BASE_CELL_OFFSET) & UInt64(0x7f))
h3_base_cell(o::H3Cell) = h3_base_cell(o.index)

# Get the 3-bit digit as position i (0-based)
@inline function h3_digit(h::UInt64, i::Int)
    0 ≤ i ≤ 15 || throw(ArgumentError("Digit index must be in 0:15"))
    shift = H3_DIGIT_BITS * (H3_NUM_DIGITS - i)
    Int((h >> shift) & H3_DIGIT_MASK)  # 0–6 when used; 7 means "unused"
end

h3_digits(x::UInt64) = h3_digit.(x, 1:h3_resolution(x))
h3_digits(o::H3Cell) = h3_digits(o.index)

@inline function h3_change_digit(x::UInt64, i::Integer, d::Integer)
    1 ≤ i ≤ H3_NUM_DIGITS || throw(ArgumentError("Digit index must be in 1:$H3_NUM_DIGITS"))
    0 ≤ d ≤ 7 || throw(ArgumentError("Digit value must be in 0:7"))
    shift = H3_DIGIT_BITS * (H3_NUM_DIGITS - i)
    mask = ~(H3_DIGIT_MASK << shift)
    (x & mask) | (UInt64(d) << shift)
end

decode(o::H3Cell) = string(h3_base_cell(o), "-", join(h3_digits(o)))

#-----------------------------------------------------------------------------# inspection
function parent(o::H3Cell)
    r = resolution(o)
    r == 0 && return nothing
    H3Cell(h3_base_cell(o.index), h3_digits(o.index)[1:end-1])
end

function children(o::H3Cell)
    child_res = Cint(resolution(o) + 1)
    size = Ref{Int64}()
    h3check(@ccall(libh3.cellToChildrenSize(o.index::UInt64, child_res::Cint, size::Ptr{Int64})::H3Error))
    out = Vector{UInt64}(undef, size[])
    h3check(@ccall(libh3.cellToChildren(o.index::UInt64, child_res::Cint, pointer(out)::Ptr{UInt64})::H3Error))
    return H3Cell.(out)
end

siblings(o::H3Cell) = resolution(o) == 0 ? nothing : filter!(!=(o), children(parent(o)))

#-----------------------------------------------------------------------------# misc
h3_n_cells(res::Integer) = 2 + 120 * 6 ^ Int(res)

h3_pentagons(r::Integer) = H3Cell.(H3_PENTAGON_BASE_CELLS, Ref(fill(0, r)))  # Faster than ccall
pentagons(o::H3Grid, r::Integer) = h3_pentagons(r)

h3_string(idx::UInt64) = string(idx, base=16)
h3_string(o::H3Cell) = h3_string(o.index)

# Faces that the cell intersects
function h3_face_numbers(idx::UInt64)
    max_faces = Ref{Cint}()
    h3check(@ccall(libh3.maxFaceCount(idx::UInt64, max_faces::Ptr{Cint})::H3Error))
    faces = Vector{Cint}(undef, max_faces[])
    h3check(@ccall(libh3.getIcosahedronFaces(idx::UInt64, pointer(faces)::Ptr{Cint})::H3Error))
    return filter!(!=(-1), faces)
end
h3_face_numbers(o::H3Cell) = h3_face_numbers(o.index)

function grid_distance(a::H3Cell, b::H3Cell)
    dist = Ref{Int64}()
    h3check(@ccall(libh3.gridDistance(a.index::UInt64, b.index::UInt64, dist::Ptr{Int64})::H3Error))
    return dist[]
end

function grid_ring(o::H3Cell, k::Integer)
    k ≥ 0 || throw(ArgumentError("k must be non-negative.  Found: $k."))
    num_cells = k == 0 ? 1 : 6k
    out = Vector{UInt64}(undef, num_cells)
    h3check(@ccall(libh3.gridRing(o.index::UInt64, Cint(k)::Cint, pointer(out)::Ptr{UInt64})::H3Error))
    filter!(!iszero, out)
    return H3Cell.(out)
end

function grid_disk(o::H3Cell, k::Integer)
    k ≥ 0 || throw(ArgumentError("k must be non-negative.  Found: $k."))
    num_cells = Cint(3 * k * (k + 1) + 1)
    out = Vector{UInt64}(undef, num_cells)
    h3check(@ccall(libh3.gridDisk(o.index::UInt64, Cint(k)::Cint, pointer(out)::Ptr{UInt64})::H3Error))
    filter!(!iszero, out)
    return H3Cell.(out)
end

function grid_path(a::H3Cell, b::H3Cell)
    size = Ref{Int64}()
    h3check(@ccall(libh3.gridPathCellsSize(a.index::UInt64, b.index::UInt64, size::Ptr{Int64})::H3Error))
    out = Vector{UInt64}(undef, size[])
    h3check(@ccall(libh3.gridPathCells(a.index::UInt64, b.index::UInt64, pointer(out)::Ptr{UInt64})::H3Error))
    filter!(!iszero, out)
    return H3Cell.(out)
end


haversine(a::H3Cell, b::H3Cell) = haversine(GI.centroid(a), GI.centroid(b))

destination(a::H3Cell, azimuth°, m) = H3Cell(destination(GI.centroid(a), azimuth°, m), resolution(a))

#-----------------------------------------------------------------------------# h3cells
"""
    h3cells(geom, res; kw...)

Convenience wrapper for `cells(H3Cell, geom, res; kw...)`.
"""
h3cells(args...; kw...) = cells(H3Cell, args...; kw...)

function cells(::Type{H3Cell}, ::GI.LineTrait, geom, res::Integer; containment = :shortest_path)
    containment in (:shortest_path, :overlap, :overlap_bbox) ||
        throw(ArgumentError("Invalid containment mode.  Expected one of `(:shortest_path, :overlap, :overlap_bbox)`.  Found: $containment."))
    coords = [LonLat(x...) for x in GI.coordinates(geom)]
    a = H3Cell(coords[1], res)
    b = H3Cell(coords[2], res)
    containment == :shortest_path && return grid_path(a, b)
    coords = [[coords..., reverse(coords)..., coords[1]]]
    h3cells(GI.Polygon(coords), res; containment)
end

function cells(::Type{H3Cell}, ::GI.PolygonTrait, geom, res::Integer; containment = :overlap)
    containment in (:center, :full, :overlap, :overlap_bbox) ||
        throw(ArgumentError("Invalid containment mode.  Expected one of `(:center, :full, :overlap, :overlap_bbox)`.  Found: $containment."))
    verts = map(GI.coordinates(geom)) do ring
        map(ring) do coord
            H3LatLng(deg2rad(coord[2]), deg2rad(coord[1]))
        end
    end
    GC.@preserve verts begin
        geo_loops = H3GeoLoop.(length.(verts), pointer.(verts))
        n = length(geo_loops)
        geo_polygon = H3GeoPolygon(geo_loops[1], n - 1, n > 1 ? pointer(geo_loops, 2) : C_NULL)
        flag_dict = Dict{Symbol, UInt32}(:center => 0, :full => 1, :overlap => 2, :overlap_bbox => 3)
        flag = flag_dict[containment]
        max_n = Ref{Int64}()
        h3check(@ccall(libh3.maxPolygonToCellsSizeExperimental(Ref(geo_polygon)::Ptr{H3GeoPolygon}, Cint(res)::Cint, flag::Int32, max_n::Ptr{Int64})::H3Error))
        out = zeros(UInt64, max_n[])
        h3check(@ccall(libh3.polygonToCellsExperimental(Ref(geo_polygon)::Ptr{H3GeoPolygon}, Cint(res)::Cint, flag::Int32, Cint(max_n[])::Cint, pointer(out)::Ptr{UInt64})::H3Error))
    end
    return H3Cell.(filter!(!iszero, unique!(out)))
end
