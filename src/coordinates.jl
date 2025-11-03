#-----------------------------------------------------------------------------# constants
"Authalic radius of the earth in meters (WGS84)"
const R = 6_371_007.180918475

"Equitorial radius (semi-major axis) of the earth (WGS84)."
const a = 6_378_137.0

"Polar radius (semi-minor axis) of the earth (WGS84)."
const b = 6_356_752.0

"Reverse flattening of the earth (WGS84)."
const rf = 298.257223563

"Flattening of the earth (WGS84)."
const f = 1 / rf

"Eccentricity of the earth (WGS84)."
const e = 0.0818191908426215

"First eccentricity squared of the earth (WGS84)."
const e2 = e ^ 2

"Second eccentricity squared of the earth (WGS84)."
const e′2 = a ^ 2 / b ^ 2 - 1.0  # second eccentricity squared

"""
    prime_vertical_radius(lat°)

Reference:  https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical
"""
prime_vertical_radius(lat) = a ^ 2 / sqrt(a ^ 2 * cosd(lat) ^ 2 + b ^ 2 * sind(lat) ^ 2)


#-----------------------------------------------------------------------------# Transformations
rotate(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
skew(φ) = [1 tan(φ); 0 1]

#-----------------------------------------------------------------------------# AbstractCoordinates
# Wrapper around a Point{N, Float64}
abstract type AbstractCoordinates{N} end
Base.getindex(o::AbstractCoordinates{N}, i::Int) where {N} = GI.coordinates(o)[i]
Base.length(o::AbstractCoordinates{N}) where {N} = N
Base.iterate(o::AbstractCoordinates{N}, i...) where {N} = iterate(o.pt, i...)

GI.isgeometry(::AbstractCoordinates) = true
GI.geomtrait(o::AbstractCoordinates) = GI.PointTrait()
GI.ncoord(::GI.PointTrait, o::AbstractCoordinates{N}) where {N} = N
GI.getcoord(::GI.PointTrait, o::AbstractCoordinates{N}, i::Int) where {N} = o.pt[i]
GI.coordinates(::GI.PointTrait, o::AbstractCoordinates{N}) where {N} = o.pt

function Base.show(io::IO, o::T) where {T <: AbstractCoordinates}
    print(io, styled"{bright_cyan:$(T.name.name)} (")
    join(io, o.pt, ", ")
    print(io, ')')
end

get_proj(A::Type, B::Type) = Proj.Transformation(proj_string(A), proj_string(B); always_xy=true)

proj_transform(o::AbstractCoordinates, T::Type) = T(get_proj(typeof(o), T)(o...)...)

#-----------------------------------------------------------------------------# Types
macro coordinates(T, N, proj_string)
    esc(quote
        """
            $($T)(pt::Point)
        """
        struct $T <: AbstractCoordinates{$N}
            pt::Point{$N, Float64}
        end
        $T(x::Real, xs::Real...) = $T(Point(x, xs...))
        $T(coord::$T) = coord
        proj_string(::Type{<: $T}) = $proj_string
    end)
end

@coordinates LonLat         2 "EPSG:4326"
@coordinates LonLatAuthalic 2 "+proj=latlon +R_a"  # FIXME...same as LonLat?

# R chosen such that Y range is [-1, 1]
@coordinates ISEA           2 "+proj=isea +R=0.6390553766196734 +orient=pole +mode=plane"
@coordinates ISEACube       2 "UNKNOWN"
@coordinates ECEF           3 "+proj=geocent +datum=WGS84 +units=m +no_defs"

#-----------------------------------------------------------------------------# LonLat
LonLat(o::AbstractCoordinates) = proj_transform(o, LonLat)
LonLat(o::ISEACube) = LonLat(ISEA(o))
LonLat(o::ECEF) = LonLat(get_proj(ECEF, LonLat)(o.pt...)[1:2]...)

#-----------------------------------------------------------------------------# LonLatAuthalic
LonLatAuthalic(o::AbstractCoordinates) = proj_transform(o, LonLatAuthalic)
LonLatAuthalic(o::ISEACube) = LonLatAuthalic(ISEA(o))
LonLatAuthalic(o::ECEF) = LonLatAuthalic(get_proj(ECEF, LonLatAuthalic)(o.pt...)[1:2]...)

#-----------------------------------------------------------------------------# ISEA
ISEA(o::AbstractCoordinates) = proj_transform(o, ISEA)
ISEA(x::ISEACube) = ISEA((CUBE_MAT_INV * x.pt)...)

const CUBE_MAT = SMatrix{2,2}(skew(deg2rad(30)) * rotate(deg2rad(60)))
const CUBE_MAT_INV = SMatrix{2,2}(inv(CUBE_MAT))

#-----------------------------------------------------------------------------# ISEACube
ISEACube(x::AbstractCoordinates) = ISEACube(ISEA(x))
ISEACube(x::ISEA) = ISEACube((CUBE_MAT * x.pt)...)

#-----------------------------------------------------------------------------# ECEF
ECEF(o::LonLat) = ECEF(get_proj(LonLat, ECEF)(o[1], o[2], 0.0)...)
ECEF(o::LonLatAuthalic) = ECEF(get_proj(LonLatAuthalic, ECEF)(o[1], o[2], 0.0)...)
ECEF(o::ISEA) = ECEF(get_proj(ISEA, ECEF)(o[1], o[2], 0.0)...)
ECEF(o::ISEACube) = ECEF(ISEA(o))


#-----------------------------------------------------------------------------# haversine
function _haversine((a_lon, a_lat), (b_lon, b_lat))
    x = sind((b_lat - a_lat) / 2) ^ 2 + cosd(a_lat) * cosd(b_lat) * sind((b_lon - a_lon) / 2) ^ 2
    return 2R * asin(min(sqrt(x), one(x)))
end

"""
    haversine(a, b)

Calculate the great-circle distance (meters) between two `AbstractCoordinates` using the Haversine formula.
"""
haversine(a::AbstractCoordinates, b::AbstractCoordinates) = _haversine(LonLat(a), LonLat(b))

#-----------------------------------------------------------------------------# destination
function _destination((lon, lat), azimuth°, m)
    δ = rad2deg(m / R)
    lat2 = asind(sind(lat) * cosd(δ) + cosd(lat) * sind(δ) * cosd(azimuth°))
    lon2 = lon + atand(sind(azimuth°) * sind(δ) * cosd(lat), cosd(δ) - sind(lat) * sind(lat2))
    (lon2, lat2)
end

"""
    destination(x, azimuth°, dist_m)

Find destination point given starting point `x::AbstractCoordinates`, `azimuth` (clockwise from North), and `dist_m` (meters).
"""
destination(o::T, azimuth°, m) where {T <: AbstractCoordinates} = T(_destination(LonLat(o), azimuth°, m))

#-----------------------------------------------------------------------------# azimuth
function _azimuth((a_lon, a_lat), (b_lon, b_lat))              # assume degrees in, degrees out
    Δlon = b_lon - a_lon

    # atan2(y, x) with degree trig
    y = cosd(b_lat) * sind(Δlon)
    x = cosd(a_lat) * sind(b_lat) - sind(a_lat) * cosd(b_lat) * cosd(Δlon)
    θ = atand(y, x)                      # [-180, 180]
    return mod(θ + 360, 360)             # [0, 360)
end

"""
    azimuth(a, b)
Calculate the azimuth (degrees clockwise from North) between two `AbstractCoordinates`.
"""
azimuth(a::AbstractCoordinates, b::AbstractCoordinates) = _azimuth(LonLat(a), LonLat(b))


#-----------------------------------------------------------------------------# coord2extent
function _coords2extent(x::LonLat; n=1000, s=1000, e=1000, w=1000)
    N = destination(x, 0, n)[2]
    E = destination(x, 90, e)[1]
    S = destination(x, 180, s)[2]
    W = destination(x, 270, w)[1]
    return Extents.Extent(X=(W, E), Y=(S, N))
end

"""
    coords2extent(x::AbstractCoordinates; n=1000, s=1000, e=1000, w=1000)
    coords2extent(x::AbstractCoordinates, dist_m::Real)

Calculate an extent around the given coordinates `x`, extending `n`, `s`, `e`, and `w` meters in each direction.
Alternatively, provide a single `dist_m` value to extend equally in all directions.
"""
coords2extent(x::AbstractCoordinates; kw...) = _coords2extent(LonLat(x); kw...)
coords2extent(x::LonLat; kw...) = _coords2extent(x; kw...)
coords2extent(x::AbstractCoordinates, dist::Real) = coords2extent(x; n=dist, s=dist, e=dist, w=dist)
