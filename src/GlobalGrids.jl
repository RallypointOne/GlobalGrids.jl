module GlobalGrids

import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import GeometryBasics as GB
import GeometryBasics: Point, Point2, Point3, Point3d, Mesh
import Proj
import Extents
import StyledStrings: @styled_str

using StaticArrays
using LinearAlgebra
using Rotations
using H3_jll: libh3

export icosahedron, LonLat, cells, CellTable,
    # H3:
    H3Grid, H3Cell, h3cells, h3bin


#-----------------------------------------------------------------------------# icosahedron
include("icosahedron.jl")

#-----------------------------------------------------------------------------# LonLat
struct LonLat{T}
    lon::T
    lat::T
end
Base.show(io::IO, o::LonLat{T}) where {T} = print(io, styled"{bright_cyan:LonLat\{$T\}} (", o.lon, ", ", o.lat, ")")
LonLat((x,y)::NTuple{2, AbstractFloat}) = LonLat(x, y)

Base.getindex(o::LonLat, i::Integer) = i == 1 ? o.lon : i == 2 ? o.lat : throw(BoundsError(o, i))
Base.length(o::LonLat) = 2
Base.iterate(o::LonLat, i = 1) = i > 2 ? nothing : (o[i], i + 1)
Base.eltype(::LonLat{T}) where {T} = T

GI.isgeometry(::LonLat) = true
GI.geomtrait(o::LonLat) = GI.PointTrait()
GI.ncoord(::GI.PointTrait, o::LonLat) = 2
GI.getcoord(::GI.PointTrait, o::LonLat, i::Int) = o[i]
GI.coordinates(::GI.PointTrait, o::LonLat) = o

"Authalic radius of the earth in meters (WGS84)"
const R = 6_371_007.180918475

function haversine((a_lon, a_lat)::LonLat, (b_lon, b_lat)::LonLat)
    x = sind((b_lat - a_lat) / 2) ^ 2 + cosd(a_lat) * cosd(b_lat) * sind((b_lon - a_lon) / 2) ^ 2
    return 2R * asin(min(sqrt(x), one(x)))
end

function destination((lon, lat)::LonLat, azimuth°::Real, m::Real)
    δ = rad2deg(m / R)
    lat2 = asind(sind(lat) * cosd(δ) + cosd(lat) * sind(δ) * cosd(azimuth°))
    lon2 = lon + atand(sind(azimuth°) * sind(δ) * cosd(lat), cosd(δ) - sind(lat) * sind(lat2))
    LonLat(lon2, lat2)
end

function azimuth((a_lon, a_lat), (b_lon, b_lat))
    Δlon = b_lon - a_lon

    # atan2(y, x) with degree trig
    y = cosd(b_lat) * sind(Δlon)
    x = cosd(a_lat) * sind(b_lat) - sind(a_lat) * cosd(b_lat) * cosd(Δlon)
    θ = atand(y, x)                      # [-180, 180]
    return mod(θ + 360, 360)             # [0, 360)
end

#-----------------------------------------------------------------------------# AbstractGrid
abstract type AbstractGrid end
Base.show(io::IO, ::T) where {T <: AbstractGrid} = print(io, styled"{bright_cyan: ⣿⣿ $(T.name.name)}")

GI.isgeometry(::AbstractGrid) = true
GI.geomtrait(o::AbstractGrid) = GI.PolyhedralSurfaceTrait()

# Interface:
# GI.crs
# GI.ncoord
# GI.npolygon
# Iterates over Cell{GridType, IndexType}

#-----------------------------------------------------------------------------# AbstractCell
abstract type AbstractCell end

# Interface for AbstractCell:
# CellType(::LonLat, ::Integer)
# grid(::Cell)::GridType
# icon(::Cell)
# decode(::Cell)
# resolution(::Cell)
# Base.getindex(::GridType, ::Integer)::Cell{GridType, IndexType}
# GI.coordinates(::Cell)
is_pentagon(::AbstractCell) = false

function Base.show(io::IO, o::T) where {T <: AbstractCell}
    print(io, styled"{bright_green: $(icon(o))}")
    print(io, styled" {bright_cyan:$(T.name.name)} {bright_yellow:$(resolution(o))}")
    print(io, styled"{bright_black: $(decode(o))}")
end

GI.isgeometry(::T) where {T <: AbstractCell} = true
GI.crs(o::T) where {T <: AbstractCell} = GI.crs(grid(o))
GI.ncoord(::GI.PolygonTrait, o::T) where {T <: AbstractCell} = GI.ncoord(grid(o))
GI.geomtrait(o::T) where {T <: AbstractCell} = GI.PolygonTrait()
GI.nhole(::GI.PolygonTrait, o::T) where {T <: AbstractCell} = 0
GI.ngeom(::GI.PolygonTrait, o::T) where {T <: AbstractCell} = 1
GI.getgeom(::GI.PolygonTrait, o::T, i::Integer) where {T <: AbstractCell} = GI.LineString(GI.coordinates(o))

# Is this faster than GO default?
# # GeometryOps optimizations for containment of point in cell
# GO._within(::GI.PointTrait, xy::LonLat, ::GI.PolygonTrait, o::T) where {T <: AbstractCell} = T(xy) == o
# GO._within(::GI.PointTrait, xy, ::GI.PolygonTrait, o::T) where {T <: AbstractCell} = T(LonLat(xy...)) == o


#-----------------------------------------------------------------------------# cells
# New methods implement: cells(::Type{T}, trait, geom, res; kw...)

cells(::Type{T}, geom, res; kw...) where {T <: AbstractCell} = cells(T, GI.geomtrait(geom), geom, res; kw...)

cells(::Type{T}, ::GI.PointTrait, geom, res) where {T <: AbstractCell} = [T(LonLat(GI.coordinates(geom)...), res)]

function cells(::Type{T}, ::GI.MultiPointTrait, geom, res) where {T <: AbstractCell}
    unique(T(LonLat(x...), res) for x in GI.coordinates(geom))
end

# !!! No fallbacks for: LineTrait, PolygonTrait

function cells(::Type{T}, ::GI.LineStringTrait, geom, res::Integer; containment = :shortest_path) where {T <: AbstractCell}
    out = T[]
    coords = [LonLat(x...) for x in GI.coordinates(geom)]
    @views for (a,b) in zip(coords[1:end-1], coords[2:end])
        line = GI.Line([a, b])
        union!(out, cells(T, line, res; containment))
    end
    return out
end

function cells(::Type{T}, ::GI.MultiPolygonTrait, geom, res::Integer; kw...) where {T <: AbstractCell}
    mapreduce(x -> cells(T, x, res; kw...), union, GI.getpolygon(geom))
end


#-----------------------------------------------------------------------------# Things without geomtraits
function cells(::Type{CT}, ::Nothing, (; X, Y)::Extents.Extent, res::Integer; kw...) where {CT <: AbstractCell}
    ls = GI.LineString([(X[1], Y[1]), (X[1], Y[2]), (X[2], Y[2]), (X[2], Y[1]), (X[1], Y[1])])
    cells(CT, GI.Polygon([ls]), res; kw...)
end

# Raster Data, e.g. (z, (x, y)) = (r, r.dims)
function cells(::Type{CT}, ::Nothing, (z, (x, y))::Tuple{Z, XY}, res::Integer; dropmissing=true) where {CT <: AbstractCell, T, Z<:AbstractMatrix{T}, XY<:Tuple}
    S = dropmissing ? Base.nonmissingtype(T) : T
    out = Dict{CT, Vector{S}}()
    for ((x, y), z) in zip(Iterators.product(x, y), z)
        if !ismissing(z) || !dropmissing
            cell = CT(LonLat(x, y), res)
            data = get!(out, cell, S[])
            push!(data, z)
        end
    end
    return out
end


#-----------------------------------------------------------------------------# h3
include("h3.jl")

#-----------------------------------------------------------------------------# plotting utils
function crosses_meridian(geom, λ0 = 180.0)
    lons = Float64[x[1] for x in GI.coordinates(geom)]
    rotated = mod.(lons .- λ0, 360)  # rotate such that λ0 is the discontinuity
    return maximum(rotated) - minimum(rotated) > 180
end










# #-----------------------------------------------------------------------------# centroid
# centroid(x) = GI.centroid(x)
# area(x) = GI.area(x)

# centroid(tri::GB.Triangle{N, T}) where {N, T} = sum(tri.points) / N

# area((a, b, c)::GB.Triangle) = norm(cross(b - a, c - a)) / 2


# #-----------------------------------------------------------------------------# Grid
# struct Grid{T}
#     ico::Mesh{3, T, GB.TriangleFace{Int}}
#     dual::Mesh{3, T, GB.NgonFace{5, Int}}
# end
# function Grid{T}() where {T}
#     # Rotate X by 30° to get pointy top (vertex at north pole)
#     ico = rotate(unit_icosahedron(T), Rotations.RotZX(deg2rad(-36), deg2rad(90)))
#     return Grid(ico, dual(ico))
# end
# Grid(T = Float64) = Grid{T}()

# Base.show(io::IO, o::Grid{T}) where {T} = print(io, styled"{bright_cyan: ⣿⣿ Grid\{$T\}}")

# #-----------------------------------------------------------------------------# face interface
# face_normals(grid::Grid{T}) where {T} = grid.dual.position
# face_index(g::Grid, pt::Point3) = findmax(x -> dot(x, pt), face_normals(g))[2]
# face_index(g::Grid, pt::AbstractCoordinates) = face_index(g, ECEF(pt).pt)
# face(g::Grid, pt)  = g.ico[face_index(g, pt)]
# faces(f, g::Grid)  = (f(tri) for tri in g.ico)
# faces(g::Grid)  = (tri for tri in g.ico)

# #-----------------------------------------------------------------------------# base cell interface
# base_normals(g::Grid) = g.ico.position
# base_cell_index(g::Grid, pt::Point) = findmax(x -> dot(x, pt), base_normals(g))[2]
# base_cell_index(g::Grid, pt::AbstractCoordinates) = base_cell_index(g, ECEF(pt).pt)
# base_cell(g::Grid, pt) = g.dual[base_cell_index(g, pt)]
# base_cells(f, g::Grid) = (f(pent) for pent in g.dual)
# base_cells(g::Grid) = (pent for pent in g.dual)

# # #-----------------------------------------------------------------------------# Cell
# # struct Cell{T}
# #     grid::Grid{T}
# #     index::UInt64
# # end

# # function encode_index(i::Integer, digits::Vector{UInt8} = UInt8[])
# #     0 < i < 21 || throw(ArgumentError("Base cell index must be in 1-20.  Found: $i."))
# #     idx = UInt64(i) << (Z7_MAX_RES * Z7_DIGIT_BITS)
# # end


# # function Base.getindex(g::Grid{T}, i::Integer) where {T}
# #     0 < i < 21 || throw(ArgumentError("Base cell index must be in 1-20.  Found: $i."))

# #     for (i, d) in enumerate(digits)
# #         idx |= (UInt64(d) << ((i - 1) * Z7_DIGIT_BITS))
# #     end
# #     for i in (length(digits)+1):Z7_MAX_RES
# #         idx |= (UInt64(Z7_PAD_VALUE) << ((i - 1) * Z7_DIGIT_BITS))
# #     end
# #     Cell{T}(g, idx)
# # end

# # function Base.show(io::IO, o::Cell{T}) where {T}
# #     print(io, styled"{bright_green: ⎔ Cell\{$T\} base=$(base(o)) res=$(resolution(o))}, {bright_black:$o}")
# # end
# # function Base.print(io::IO, o::Cell{T}) where {T}
# #     b = base(o)
# #     d = digits(o)
# #     filter!(!=(0x07), d)
# #     print(io, lpad(Int(b), 2, '0'), '-')
# #     join(io, Int.(d))
# # end

# # function Base.digits!(x::AbstractVector{UInt8}, o::Cell)
# #     mask = (UInt64(1) << Z7_DIGIT_BITS) - 1  # 0b111
# #         for j in 0:(Z7_MAX_RES-1)
# #             d = UInt8((o.index >> (j * Z7_DIGIT_BITS)) & mask)
# #             d === Z7_PAD_VALUE && break
# #             push!(x, d)
# #         end
# #     return x
# # end
# # Base.digits(o::Cell) = digits!(UInt8[], o)

# # base(o::Cell) = o.index >> (Z7_MAX_RES * Z7_DIGIT_BITS)

# # resolution(o::Cell) = length(digits(o))




# # #-----------------------------------------------------------------------------# Z7 index
# # const Z7_BASE_BITS  = 4
# # const Z7_DIGIT_BITS = 3
# # const Z7_MAX_RES    = 20
# # const Z7_PAD_VALUE  = 0x07

# # # (q, r) axial coordinates for the 7 child cells of a parent
# # const Z7_OFFSETS = (
# #     ( 0,  0),  # digit 0 (center)
# #     ( 1,  0),  # 1
# #     ( 0,  1),  # 2
# #     (-1,  1),  # 3
# #     (-1,  0),  # 4
# #     ( 0, -1),  # 5
# #     ( 1, -1),  # 6
# # )









# # function Index(base, digits)
# #         -1 < base < 12 || throw(ArgumentError("base must be in 0-11.  Found: $base"))
# #     all(x -> (-1 < x < 7), digits) || throw(ArgumentError("digits must all be in 0-6.  Found: $digits"))
# #     idx = UInt64(base - 1) << (Z7_MAX_RES * Z7_DIGIT_BITS)
# #     for (i, d) in enumerate(digits)
# #         idx |= (UInt64(d) << ((i - 1) * Z7_DIGIT_BITS))
# #     end
# #     for i in (length(digits)+1):Z7_MAX_RES
# #         idx |= (UInt64(Z7_PAD_VALUE) << ((i - 1) * Z7_DIGIT_BITS))
# #     end
# #     return Index(idx)
# # end





# # # Every location on the grid can be represented by (face, q, r)

# # function face_xy(z::Index)
# #     face = base(z)
# #     q, r = 0.0, 0.0
# #     scale = 1.0
# #     for d in digits(z)
# #         dq, dr = Z7_OFFSETS[d + 1]
# #         scale /= sqrt(7)
# #         q += dq * scale
# #         r += dr * scale
# #     end
# #     (face, q, r)
# # end


# # Base.getindex(o::Grid{T}, coord::AbstractCoordinates) where {T} = getindex(o, ECEF(coord))

# # function Base.getindex(o::Grid{T}, coord::ECEF) where {T}
# #     v = coord.pt
# #     f = base(o.ico, v)
# #     v0, v1, v2 = o.ico[f]
# #     p = cross(cross(v0, v), cross(v1, v2))
# #     h = sqrt((1 - dot(v0, v)) / (1 - dot(v0, p)))
# #     β2 = h * area(GB.Triangle(v0, v1, p)) / area(GB.Triangle(v0, v1, v2))
# #     β0 = 1 - h
# #     β1 = h - β2
# #     β0, β1, β2
# # end

# # #-----------------------------------------------------------------------------# Cell
# # struct Cell{T}
# #     grid::Grid{T}
# #     resolution::Int
# #     base::Int
# #     digits::Vector{Int}
# #     center::Point3{T}
# # end
# # function Base.show(io::IO, o::Cell{T}) where {T}
# #     content = [
# #         styled" {bright_yellow:res=$(o.resolution)}",
# #         styled" {bright_green:base=$(o.base)}",
# #         styled" {bright_cyan:$(o.digits)}",
# #         styled" {bright_black:$(o.center)}"
# #     ]
# #     print(io, "Cell{$T}", content...)
# # end

# # # center of cell `z`
# # function centroid(grid::Grid{T}, z::Index) where {T}
# #     f, q, r = face_xy(z)
# #     pt = ISEA{T}(Point2{T}(q, r))
# #     pt
# # end

# # function cell(grid::Grid{T}, pt::Point3{T}, res::Integer) where {T}
# #     i = base(grid.ico, pt)
# # end

# # #-----------------------------------------------------------------------------# getindex methods for Grid and Cell
# # # grid[::Integer]
# # function Base.getindex(grid::Grid{T}, i::Integer) where {T}
# #     0 < i < 21 || throw(ArgumentError("Base cell index must be in 1-20.  Found: $i"))
# #     return Cell{T}(grid, 0, i, Int[], Point3{T}(grid.ico.face_normals[i]))
# # end

# # function Base.getindex(o::Cell{T}, i::Integer) where {T}
# #     0 < i < 8 || throw(ArgumentError("Child digit must be in 0-6.  Found: $i"))
# #     # TODO: compute actual center point
# #     center = Point3{T}((0,0,0))
# #     return Cell{T}(o.grid, o.resolution + 1, o.base, vcat(o.digits, i), center)
# # end

end # module GlobalGrids
