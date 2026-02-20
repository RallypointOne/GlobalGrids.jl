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

export icosahedron, LonLat, cells,
    # H3:
    H3Grid, H3Cell, h3cells


#-----------------------------------------------------------------------------# icosahedron
include("icosahedron.jl")

#-----------------------------------------------------------------------------# LonLat
"""
    LonLat(lon, lat)
    LonLat((lon, lat))

A point on the Earth's surface in longitude/latitude (degrees).  Implements the
GeoInterface `PointTrait` so it works with plotting and spatial libraries.

### Examples

    LonLat(-75.0, 54.0)
"""
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

"Authalic radius of the earth in meters (WGS84)."
const EARTH_RADIUS_AUTHALIC = 6_371_007.180918475

"""
    haversine(a::LonLat, b::LonLat)

Great-circle distance (meters) between two points using the Haversine formula.
"""
function haversine((a_lon, a_lat)::LonLat, (b_lon, b_lat)::LonLat)
    x = sind((b_lat - a_lat) / 2) ^ 2 + cosd(a_lat) * cosd(b_lat) * sind((b_lon - a_lon) / 2) ^ 2
    return 2EARTH_RADIUS_AUTHALIC * asin(min(sqrt(x), one(x)))
end

"""
    destination(origin::LonLat, azimuth°::Real, m::Real)

Destination point given a starting `origin`, `azimuth°` (degrees clockwise from
North), and distance `m` (meters).
"""
function destination((lon, lat)::LonLat, azimuth°::Real, m::Real)
    δ = rad2deg(m / EARTH_RADIUS_AUTHALIC)
    lat2 = asind(sind(lat) * cosd(δ) + cosd(lat) * sind(δ) * cosd(azimuth°))
    lon2 = lon + atand(sind(azimuth°) * sind(δ) * cosd(lat), cosd(δ) - sind(lat) * sind(lat2))
    LonLat(lon2, lat2)
end

"""
    azimuth(a, b)

Initial bearing (degrees clockwise from North) from point `a` to point `b`.
"""
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

#-----------------------------------------------------------------------------# AbstractCell
abstract type AbstractCell end

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

#-----------------------------------------------------------------------------# cells
"""
    cells(::Type{T}, geom, res; kw...) where {T <: AbstractCell}

Return the grid cells of type `T` at resolution `res` that cover `geom`.

Dispatches on the GeoInterface geometry trait of `geom` and supports points,
multipoints, lines, linestrings, polygons, multipolygons, and `Extents.Extent`.
"""
cells(::Type{T}, geom, res; kw...) where {T <: AbstractCell} = cells(T, GI.geomtrait(geom), geom, res; kw...)

cells(::Type{T}, ::GI.PointTrait, geom, res) where {T <: AbstractCell} = [T(LonLat(GI.coordinates(geom)...), res)]

function cells(::Type{T}, ::GI.MultiPointTrait, geom, res) where {T <: AbstractCell}
    unique(T(LonLat(x...), res) for x in GI.coordinates(geom))
end

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

end # module GlobalGrids
