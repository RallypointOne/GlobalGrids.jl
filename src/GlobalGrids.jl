module GlobalGrids

import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import GeometryBasics as GB
import GeometryBasics: Point, Point2, Point3
import Proj
import StyledStrings: @styled_str
import StaticArrays as SA

using LinearAlgebra

export UnitIcosahedron, Grid, Index, Cell, base, resolution

#-----------------------------------------------------------------------------# constants
include("GeoConstants.jl")
using .GeoConstants: R

#-----------------------------------------------------------------------------# Coordinate Wrappers
abstract type CoordinateWrapper end
Base.getindex(o::CoordinateWrapper, i::Int) = o.pt[i]
Base.iterate(o::CoordinateWrapper, x...) = iterate(o.pt, x...)

# Latitude is Authalic!!!!
struct LonLatRad{T} <: CoordinateWrapper
    pt::Point2{T}
end
struct LonLatDeg{T} <: CoordinateWrapper
    pt::Point2{T}
end
struct ISEA{T} <: CoordinateWrapper
    pt::Point2{T}
end
struct ECEF{T} <: CoordinateWrapper
    pt::Point3{T}
end

LonLatRad(o::LonLatRad) = o
LonLatRad(o::LonLatDeg) = LonLatRad(deg2rad.(o.pt))
LonLatRad(o::ISEA) = LonLatRad(LonLatDeg(o))
LonLatRad((; pt)::ECEF) = (x = atan(pt[2], pt[1]); hyp = norm(pt); y = atan(pt[3], hyp); LonLatRad(Point2(x, y)))

LonLatDeg(o::LonLatRad) = LonLatDeg(rad2deg.(o.pt))
LonLatDeg(o::LonLatDeg) = o
function LonLatDeg(o::ISEA)
    proj = Proj.Transformation("+proj=isea +R=$R +orient=isea +mode=plane", "EPSG:4326", ; always_xy=true)
    x, y = proj(o.pt...)
    return LonLatDeg(Point2(x, y))
end
LonLatDeg(o::ECEF) = LonLatDeg(LonLatRad(o))

ISEA(o::LonLatRad) = ISEA(LonLatDeg(o))
function ISEA(o::LonLatDeg)
    proj = Proj.Transformation("EPSG:4326", "+proj=isea +R=$R +orient=isea +mode=plane"; always_xy=true)
    x, y = proj(o.pt...)
    return ISEA(Point2(x, y))
end
ISEA(o::ISEA) = o
ISEA(o::ECEF) = ISEA(LonLatDeg(o))

function ECEF((lon, lat)::LonLatRad)
    x = sin(lat) * cos(lon)
    y = sin(lat) * sin(lon)
    z = cos(lat)
    ECEF(Point3(x, y, z))
end
ECEF(o::LonLatDeg) = ECEF(LonLatRad(o))
ECEF(o::ISEA) = ECEF(LonLatDeg(o))
ECEF(o::ECEF) = o

#-----------------------------------------------------------------------------# centroid
centroid(x) = GI.centroid(x)

function centroid(tri::GB.Triangle{3, T}) where {T}
    pts = tri.points
    return GB.Point{3,T}(
        sum(p[1] for p in pts) / 3,
        sum(p[2] for p in pts) / 3,
        sum(p[3] for p in pts) / 3
    )
end
function centroid(tri::GB.Triangle{2, T}) where {T}
    pts = tri.points
    return GB.Point{2,T}(
        sum(p[1] for p in pts) / 3,
        sum(p[2] for p in pts) / 3,
    )
end

#-----------------------------------------------------------------------------# UnitUnitIcosahedron
# Unit UnitIcosahedron in Cartesian coordinates: 12 vertices, 20 triangular faces
struct UnitIcosahedron{T}
    mesh::GB.Mesh{3, T, GB.TriangleFace{Int}}
    face_normals::Vector{GB.Vec{3, T}}
    function UnitIcosahedron{T}() where {T}
        φ = (1 + sqrt(5)) / 2
        verts = map(normalize, GB.Point{3,T}[
            (-1,  φ, 0),  (1,  φ, 0),  (-1, -φ, 0),  (1, -φ, 0),
            (0, -1,  φ), (0,  1,  φ), (0, -1, -φ), (0,  1, -φ),
            (φ,  0, -1), (φ,  0,  1), (-φ, 0, -1), (-φ, 0,  1)
        ])
        faces = GB.TriangleFace.([(1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12), (2,6,10), (6,12,5),
            (12,11,3), (11,8,7), (8,2,9), (4,10,5), (4,5,3), (4,3,7), (4,7,9),  (4,9,10), (5,10,6),
            (3,5,12), (7,3,11), (9,7,8), (10,9,2)
        ])
        return new{T}(GB.Mesh(verts, faces), values(GB.face_normals(verts, faces)))
    end
end
UnitIcosahedron() = UnitIcosahedron{Float64}()

# Iterators over faces, optionally applying a function
faces(ico::UnitIcosahedron{T}) where {T} = (tri for tri in ico)
faces(f, ico::UnitIcosahedron{T}) where {T} = (f(tri) for tri in ico)

face_centers(ico::UnitIcosahedron{T}, S=ECEF) where {T} = faces(ico) do tri
    S(ECEF(centroid(tri)))
end


# function to_lonlat(ico::UnitIcosahedron{T}) where {T}
#     coords = GO.reproject(ico.mesh.position, ico.proj, LONLAT_PROJ)
#     pts = [GB.Point2{T}(x[1:2]) for x in coords]
#     mesh = GB.Mesh(pts, ico.mesh.faces)
#     # UnitIcosahedron{T, typeof(LONLAT_PROJ)}(mesh, values(GB.face_normals(pts, ico.mesh.faces)), LONLAT_PROJ)
# end


# # GeometryOps
# function reproject(ico::UnitIcosahedron{T}, proj::P) where {T, P}
#     coords = GB.Point3{T}.(GO.reproject(ico.mesh.position, ico.proj, proj))
#     mesh = GB.Mesh(coords, ico.mesh.faces)
#     UnitIcosahedron{T, P}(mesh, values(GB.face_normals(coords, ico.mesh.faces)), proj)
# end

# # GeoInterface methods
# GI.isgeometry(::UnitIcosahedron) = true
# GI.trait(o::UnitIcosahedron) = GI.PolyhedralSurfaceTrait()
# GI.ncoord(::GI.PolyhedralSurfaceTrait, o::UnitIcosahedron) = 3
# GI.ngeom(::GI.PolyhedralSurfaceTrait, o::UnitIcosahedron) = 20
# GI.getgeom(::GI.PolyhedralSurfaceTrait, o::UnitIcosahedron, i::Int) = o.mesh[i]
# GI.crs(o::UnitIcosahedron) = o.proj

# # GeometryBasics methods
# GB.coordinates(ico::UnitIcosahedron{T}) where {T} = GB.coordinates(ico.mesh)  # vertices
# GB.face_normals(ico::UnitIcosahedron{T}) where {T} = ico.face_normals

# Iteration
Base.getindex(ico::UnitIcosahedron{T}, i::Int) where {T} = ico.mesh[i]
Base.iterate(ico::UnitIcosahedron{T}, i=1) where {T} = i > 20 ? nothing : (ico[i], i + 1)
Base.length(ico::UnitIcosahedron{T}) where {T} = 20

# # Circumscribed sphere
GB.Sphere(ico::UnitIcosahedron{T}) where {T} = GB.Sphere(GB.Point{3,T}(0, 0, 0), T(1))

Base.show(io::IO, o::UnitIcosahedron{T}) where {T} = print(io, "UnitIcosahedron{$T}()")

# Get the face index (base pentagon cell)
base(ico::UnitIcosahedron{T}, pt::GB.Point{3, T}) where {T} = findmax(x -> dot(x, pt), values(ico.face_normals))[2]



# function to_lonlat(tri::GB.Triangle{3, T}) where {T}
#     proj = Proj.Transformation("EPSG:4326", ISEA_PROJ; always_xy=true)
#     coords = map(tri.points) do (x, y, z)
#         lon = atan(y, x)
#         hyp = sqrt(x^2 + y^2)
#         lat = atan(z, hyp)
#         GB.Point{2,T}(proj(lon, lat)...)
#     end
#     GB.Triangle(coords...)
# end

# function to_isea(tri::GB.Triangle{3, T}) where {T}
#     proj = Proj.Transformation("EPSG:4326", ISEA_PROJ; always_xy=true)
#     coords = map(tri.points) do (x, y, z)
#         lon = atan(y, x)
#         hyp = sqrt(x^2 + y^2)
#         lat = atan(z, hyp)
#         GB.Point{2,T}(proj(lon, lat)...)
#     end
#     GB.Triangle(coords...)
# end

#-----------------------------------------------------------------------------# Grid
struct Grid{T}
    ico::UnitIcosahedron{T}
    Grid{T}() where {T} = new{T}(UnitIcosahedron{T}())
end
Grid() = Grid{Float64}()

#-----------------------------------------------------------------------------# Index
const Z7_BASE_BITS  = 4
const Z7_DIGIT_BITS = 3
const Z7_MAX_RES    = 20
const Z7_PAD_VALUE  = 0x07

# (q, r) axial coordinates for the 7 child cells of a parent
const Z7_OFFSETS = (
    ( 0,  0),  # digit 0 (center)
    ( 1,  0),  # 1
    ( 0,  1),  # 2
    (-1,  1),  # 3
    (-1,  0),  # 4
    ( 0, -1),  # 5
    ( 1, -1),  # 6
)

struct Index
    value::UInt64
end
function Index(base, digits)
        -1 < base < 12 || throw(ArgumentError("base must be in 0-11.  Found: $base"))
    all(x -> (-1 < x < 7), digits) || throw(ArgumentError("digits must all be in 0-6.  Found: $digits"))
    idx = UInt64(base - 1) << (Z7_MAX_RES * Z7_DIGIT_BITS)
    for (i, d) in enumerate(digits)
        idx |= (UInt64(d) << ((i - 1) * Z7_DIGIT_BITS))
    end
    for i in (length(digits)+1):Z7_MAX_RES
        idx |= (UInt64(Z7_PAD_VALUE) << ((i - 1) * Z7_DIGIT_BITS))
    end
    return Index(idx)
end

function Base.print(io::IO, i::Index)
    b = base(i)
    d = digits(i)
    filter!(!=(0x07), d)
    print(io, lpad(Int(b), 2, '0'), '_')
    join(io, Int.(d))
end

Base.show(io::IO, i::Index) = print(io, "Index: ", i)

base(i::Index) = i.value >> (Z7_MAX_RES * Z7_DIGIT_BITS)

Base.digits(i::Index) = digits!(UInt8[], i)

function Base.digits!(x::AbstractVector{UInt8}, i::Index)
    mask = (UInt64(1) << Z7_DIGIT_BITS) - 1  # 0b111
    for j in 0:(Z7_MAX_RES-1)
        d = UInt8((i.value >> (j * Z7_DIGIT_BITS)) & mask)
        d === Z7_PAD_VALUE && break
        push!(x, d)
    end
    return x
end
resolution(i::Index) = length(digits(i))

Base.getindex(o::Index, i::Integer) = Index(base(o), vcat(digits(o), i))

# Every location on the grid can be represented by (face, q, r)

function face_xy(z::Index)
    face = base(z)
    q, r = 0.0, 0.0
    scale = 1.0
    for d in digits(z)
        dq, dr = Z7_OFFSETS[d + 1]
        scale /= sqrt(7)
        q += dq * scale
        r += dr * scale
    end
    (face, q, r)
end

#-----------------------------------------------------------------------------# Cell
struct Cell{T}
    grid::Grid{T}
    resolution::Int
    base::Int
    digits::Vector{Int}
    center::GB.Point3{T}
end
function Base.show(io::IO, o::Cell{T}) where {T}
    content = [
        styled" {bright_yellow:$(o.resolution)}",
        styled" {bright_green:$(o.base)}",
        styled" {bright_cyan:$(o.digits)}",
        styled" {bright_black:$(o.center)}"
    ]
    print(io, "Cell{$T}", content...)
end

# center of cell `z`
function centroid(grid::Grid{T}, z::Index) where {T}
    f, q, r = face_xy(z)
    pt = ISEA{T}(Point2{T}(q, r))
    pt
end

function cell(grid::Grid{T}, pt::GB.Point3{T}, res::Integer) where {T}
    i = base(grid.ico, pt)
end

#-----------------------------------------------------------------------------# getindex methods for Grid and Cell
# grid[::Integer]
function Base.getindex(grid::Grid{T}, i::Integer) where {T}
    0 < i < 21 || throw(ArgumentError("Base cell index must be in 1-20.  Found: $i"))
    return Cell{T}(grid, 0, i, Int[], GB.Point3{T}(grid.ico.face_normals[i]))
end

function Base.getindex(o::Cell{T}, i::Integer) where {T}
    0 < i < 8 || throw(ArgumentError("Child digit must be in 0-6.  Found: $i"))
    # TODO: compute actual center point
    center = GB.Point3{T}((0,0,0))
    return Cell{T}(o.grid, o.resolution + 1, o.base, vcat(o.digits, i), center)
end

end # module GlobalGrids
