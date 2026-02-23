#-----------------------------------------------------------------------------# Constants
const IGEO7_NUM_DIGITS  = 20
const IGEO7_DIGIT_BITS  = 3
const IGEO7_DIGIT_MASK  = UInt64(0x7)
const IGEO7_BASE_OFFSET = 60  # bits [63:60]
const IGEO7_MAX_RES     = 20

# Aperture-7 rotation angle per level (child grid rotates clockwise)
const IGEO7_THETA7 = atan(sqrt(3) / 5)  # ≈ 19.1066°

# Digit offsets in axial (q, r) coordinates
const IGEO7_DIGIT_OFFSETS = SVector{7,SVector{2,Int}}(
    SVector(0,  0),  # digit 0: center
    SVector(1,  0),  # digit 1
    SVector(0,  1),  # digit 2
    SVector(-1, 1),  # digit 3
    SVector(-1, 0),  # digit 4
    SVector(0, -1),  # digit 5
    SVector(1, -1),  # digit 6
)

#-----------------------------------------------------------------------------# ISEA Icosahedron Geometry
# 12 vertices of icosahedron(:isea) as unit vectors (1-indexed, matching GeometryBasics)
const IGEO7_VERTICES_3D = SVector{12,SVector{3,Float64}}(
    SVector(0.5181369092970741,  0.10306383925506017,  0.8490653615959628),
    SVector(0.9966866752635933, -0.0015397574699671436, -0.08132220175904135),
    SVector(-0.9966866752635933, 0.0015397574699671436, 0.08132220175904135),
    SVector(-0.5181369092970741, -0.10306383925506017, -0.8490653615959628),
    SVector(-0.4613784348709131, -0.8762800214612103, -0.13879216052791962),
    SVector(0.4748340273475299,  -0.8135346882413987,  0.3356992068815409),
    SVector(-0.4748340273475299, 0.8135346882413987,  -0.3356992068815409),
    SVector(0.4613784348709131,  0.8762800214612103,   0.13879216052791962),
    SVector(0.3829968865764098,  0.4375553752195837,  -0.8135469737447408),
    SVector(0.39131290006572633, -0.606807550066165,   -0.6918517264738518),
    SVector(-0.39131290006572633, 0.606807550066165,    0.6918517264738518),
    SVector(-0.3829968865764098, -0.4375553752195837,   0.8135469737447408),
)

# Vertex neighbors (1-indexed, sorted)
const IGEO7_VERTEX_NEIGHBORS = SVector{12,SVector{5,Int}}(
    SVector(2, 6, 8, 11, 12),   # vertex 1
    SVector(1, 6, 8, 9, 10),    # vertex 2
    SVector(4, 5, 7, 11, 12),   # vertex 3
    SVector(3, 5, 7, 9, 10),    # vertex 4
    SVector(3, 4, 6, 10, 12),   # vertex 5
    SVector(1, 2, 5, 10, 12),   # vertex 6
    SVector(3, 4, 8, 9, 11),    # vertex 7
    SVector(1, 2, 7, 9, 11),    # vertex 8
    SVector(2, 4, 7, 8, 10),    # vertex 9
    SVector(2, 4, 5, 6, 9),     # vertex 10
    SVector(1, 3, 7, 8, 12),    # vertex 11
    SVector(1, 3, 5, 6, 11),    # vertex 12
)

# Hex size parameter at resolution 1 (Lambert equal-area projection, unit sphere)
# Derived from: hex_area = 4π / n_cells, d² = 2*hex_area/√3, s = d/√3
const IGEO7_S1 = sqrt(2 * 4π / (72 * sqrt(3))) / sqrt(3)  # ≈ 0.2592

#-----------------------------------------------------------------------------# Local Basis
# Precompute local tangent-plane basis (u, v) for each base vertex.
# u is oriented toward the angular midpoint of the first two neighbors.
function _compute_local_bases()
    bases_u = Vector{SVector{3,Float64}}(undef, 12)
    bases_v = Vector{SVector{3,Float64}}(undef, 12)
    for i in 1:12
        n = IGEO7_VERTICES_3D[i]
        # Use direction toward first neighbor as reference
        nb = IGEO7_VERTICES_3D[IGEO7_VERTEX_NEIGHBORS[i][1]]
        t = nb - dot(nb, n) * n  # project onto tangent plane
        u = normalize(t)
        v = cross(n, u)
        bases_u[i] = u
        bases_v[i] = v
    end
    return SVector{12}(bases_u), SVector{12}(bases_v)
end

const IGEO7_BASIS_U, IGEO7_BASIS_V = _compute_local_bases()

#-----------------------------------------------------------------------------# IGEO7Grid
"""
    IGEO7Grid()

An aperture-7 hexagonal discrete global grid based on the ISEA icosahedron.
Uses 12 pentagonal base cells at icosahedron vertices with Z7 hierarchical
indexing.  CRS is EPSG:4326 (WGS84 longitude/latitude).

Cell count at resolution `r`: `10 * 7^r + 2`.
"""
struct IGEO7Grid <: AbstractGrid end

GI.crs(::IGEO7Grid) = GFT.EPSG(4326)
GI.ncoord(::IGEO7Grid) = 2

#-----------------------------------------------------------------------------# IGEO7Cell
"""
    IGEO7Cell(index::UInt64)
    IGEO7Cell(base::Integer, digits::AbstractVector{<:Integer})
    IGEO7Cell(coord::LonLat, res::Integer)
    IGEO7Cell(str::AbstractString)

A single cell in the IGEO7 grid.  Wraps a `UInt64` index encoding a 4-bit
base cell (0–11) and up to 20 three-bit digits (0–6, padding = 7).

### Examples

    IGEO7Cell(LonLat(-75.0, 54.0), 3)
    IGEO7Cell(0, [1, 3, 0])
"""
struct IGEO7Cell <: AbstractCell
    index::UInt64
end

#-----------------------------------------------------------------------------# Bit Manipulation
"""Extract 0-indexed base cell (0–11) from raw index."""
igeo7_base_cell(idx::UInt64) = Int((idx >> IGEO7_BASE_OFFSET) & 0xF)

"""Get 3-bit digit at 1-based position `i` (1 = most significant)."""
@inline function igeo7_digit(idx::UInt64, i::Int)
    shift = IGEO7_DIGIT_BITS * (IGEO7_NUM_DIGITS - i)
    Int((idx >> shift) & IGEO7_DIGIT_MASK)
end

"""Resolution = number of leading non-padding (non-7) digits."""
function igeo7_resolution(idx::UInt64)
    for i in 1:IGEO7_NUM_DIGITS
        igeo7_digit(idx, i) == 7 && return i - 1
    end
    return IGEO7_NUM_DIGITS
end

"""Return active digits as a vector."""
function igeo7_digits(idx::UInt64)
    r = igeo7_resolution(idx)
    [igeo7_digit(idx, i) for i in 1:r]
end

"""Hex string representation of raw index."""
igeo7_string(idx::UInt64) = string(idx, base=16)
igeo7_string(o::IGEO7Cell) = igeo7_string(o.index)

"""Number of cells at resolution `r`."""
igeo7_n_cells(r::Integer) = 10 * big(7)^r + 2

#-----------------------------------------------------------------------------# Constructors
function IGEO7Cell(base::Integer, digits::AbstractVector{<:Integer})
    0 <= base <= 11 || throw(ArgumentError("base must be in 0:11, got $base"))
    length(digits) <= IGEO7_NUM_DIGITS || throw(ArgumentError("digits length must be <= $IGEO7_NUM_DIGITS"))
    all(d -> 0 <= d <= 6, digits) || throw(ArgumentError("each digit must be in 0:6"))
    # Pentagon cells cannot have digit 6
    pent = true
    for d in digits
        if pent && d == 6
            throw(ArgumentError("pentagon cell cannot have digit 6"))
        end
        pent = (d == 0) && pent  # only center child of pentagon is pentagon
    end
    idx = UInt64(base) << IGEO7_BASE_OFFSET
    for i in 1:IGEO7_NUM_DIGITS
        d = i <= length(digits) ? UInt64(digits[i]) : IGEO7_DIGIT_MASK  # pad with 7
        shift = IGEO7_DIGIT_BITS * (IGEO7_NUM_DIGITS - i)
        idx |= (d & IGEO7_DIGIT_MASK) << shift
    end
    return IGEO7Cell(idx)
end

function IGEO7Cell(str::AbstractString)
    idx = parse(UInt64, str; base=16)
    return IGEO7Cell(idx)
end

#-----------------------------------------------------------------------------# Lambert Azimuthal Equal-Area Projection
"""Project 3D unit vector `p` onto tangent plane at base vertex `bi` (1-indexed).
Returns (x, y) in the local coordinate system."""
function _lambert_forward(p::SVector{3,Float64}, bi::Int)
    n = IGEO7_VERTICES_3D[bi]
    cos_a = clamp(dot(p, n), -1.0, 1.0)
    cos_a <= -1.0 + 1e-14 && return (0.0, 0.0)  # antipodal
    k = sqrt(2.0 / (1.0 + cos_a))
    x = k * dot(p, IGEO7_BASIS_U[bi])
    y = k * dot(p, IGEO7_BASIS_V[bi])
    return (x, y)
end

"""Inverse Lambert: (x, y) in tangent plane at base vertex `bi` → 3D unit vector."""
function _lambert_inverse(x::Float64, y::Float64, bi::Int)
    rho2 = x * x + y * y
    rho2 < 1e-30 && return IGEO7_VERTICES_3D[bi]
    cos_a = 1.0 - rho2 / 2.0
    scale = sqrt(max(0.0, 1.0 - rho2 / 4.0))
    p = cos_a * IGEO7_VERTICES_3D[bi] + scale * (x * IGEO7_BASIS_U[bi] + y * IGEO7_BASIS_V[bi])
    return normalize(p)
end

"""Convert 3D unit vector to LonLat."""
function _xyz_to_lonlat(p::SVector{3,Float64})
    lat = asind(clamp(p[3], -1.0, 1.0))
    lon = atand(p[2], p[1])
    return LonLat(lon, lat)
end

"""Convert LonLat to 3D unit vector."""
function _lonlat_to_xyz(ll::LonLat)
    clat = cosd(ll.lat)
    SVector(clat * cosd(ll.lon), clat * sind(ll.lon), sind(ll.lat))
end

#-----------------------------------------------------------------------------# Hex Math
"""Axial (q, r) → Cartesian (x, y) for pointy-top hex with size parameter `s`."""
function _axial_to_cart(q, r, s)
    x = s * (sqrt(3) * q + sqrt(3) / 2 * r)
    y = s * (1.5 * r)
    return (x, y)
end

"""Cartesian (x, y) → fractional axial (q, r) for pointy-top hex with size parameter `s`."""
function _cart_to_axial(x, y, s)
    r = y / (1.5 * s)
    q = (x / s - sqrt(3) / 2 * r) / sqrt(3)
    return (q, r)
end

"""Round fractional axial (q, r) to nearest hex center."""
function _hex_round(q::Float64, r::Float64)
    s = -q - r  # cube z
    rq = round(q)
    rr = round(r)
    rs = round(s)
    dq = abs(rq - q)
    dr = abs(rr - r)
    ds = abs(rs - s)
    if dq > dr && dq > ds
        rq = -rr - rs
    elseif dr > ds
        rr = -rq - rs
    end
    return (Int(rq), Int(rr))
end

"""Apply 2D rotation by angle `theta` (radians) to point (x, y)."""
function _rotate2d(x, y, theta)
    c, s = cos(theta), sin(theta)
    return (c * x - s * y, s * x + c * y)
end

#-----------------------------------------------------------------------------# Forward Transform (LonLat → IGEO7Cell)
function IGEO7Cell(coord::LonLat, res::Integer)
    0 <= res <= IGEO7_MAX_RES || throw(ArgumentError("resolution must be in 0:$IGEO7_MAX_RES, got $res"))
    p = _lonlat_to_xyz(coord)

    # Find nearest base vertex (0-indexed for the cell, 1-indexed for arrays)
    best_dot = -Inf
    base1 = 1  # 1-indexed
    for i in 1:12
        d = dot(p, IGEO7_VERTICES_3D[i])
        if d > best_dot
            best_dot = d
            base1 = i
        end
    end

    res == 0 && return IGEO7Cell(base1 - 1, Int[])

    # Project onto Lambert tangent plane at base vertex
    px, py = _lambert_forward(p, base1)

    # Top-down recursive digit extraction:
    # At each level, find the nearest child center among the valid digit offsets.
    digits = Vector{Int}(undef, res)
    center_x, center_y = 0.0, 0.0  # current cell center in Lambert plane
    is_pent = true  # base cells are always pentagons

    for level in 1:res
        s_level = IGEO7_S1 / sqrt(7.0)^(level - 1)
        rotation = -level * IGEO7_THETA7

        n_children = is_pent ? 6 : 7
        best_d2 = Inf
        best_digit = 0

        for d in 0:(n_children - 1)
            off = IGEO7_DIGIT_OFFSETS[d + 1]
            # Convert offset to cartesian, apply rotation
            ox, oy = _axial_to_cart(off[1], off[2], s_level)
            ox, oy = _rotate2d(ox, oy, rotation)
            cx = center_x + ox
            cy = center_y + oy
            dist2 = (px - cx)^2 + (py - cy)^2
            if dist2 < best_d2
                best_d2 = dist2
                best_digit = d
            end
        end

        # Update center and pentagon status
        off = IGEO7_DIGIT_OFFSETS[best_digit + 1]
        ox, oy = _axial_to_cart(off[1], off[2], s_level)
        ox, oy = _rotate2d(ox, oy, rotation)
        center_x += ox
        center_y += oy
        digits[level] = best_digit
        is_pent = is_pent && (best_digit == 0)
    end

    return IGEO7Cell(base1 - 1, digits)
end

#-----------------------------------------------------------------------------# Inverse Transform (IGEO7Cell → center LonLat)
"""Compute the Lambert-plane center (x, y) for a cell at its base vertex."""
function _igeo7_center_lambert(idx::UInt64)
    base1 = igeo7_base_cell(idx) + 1
    res = igeo7_resolution(idx)
    cx, cy = 0.0, 0.0
    for level in 1:res
        d = igeo7_digit(idx, level)
        s_level = IGEO7_S1 / sqrt(7.0)^(level - 1)
        rotation = -level * IGEO7_THETA7
        off = IGEO7_DIGIT_OFFSETS[d + 1]
        ox, oy = _axial_to_cart(off[1], off[2], s_level)
        ox, oy = _rotate2d(ox, oy, rotation)
        cx += ox
        cy += oy
    end
    return base1, cx, cy
end

function _igeo7_centroid(idx::UInt64)
    base1, cx, cy = _igeo7_center_lambert(idx)
    p = _lambert_inverse(cx, cy, base1)
    return _xyz_to_lonlat(p)
end

#-----------------------------------------------------------------------------# Cell Boundary
function _igeo7_boundary(idx::UInt64)
    base1, cx, cy = _igeo7_center_lambert(idx)
    res = igeo7_resolution(idx)
    pent = _igeo7_is_pentagon(idx)
    n_verts = pent ? 5 : 6

    # Hex vertex radius = circumradius = size parameter at this resolution
    s_res = res == 0 ? IGEO7_S1 : IGEO7_S1 / sqrt(7.0)^(res - 1)
    # For pointy-top hex, vertices are at distance s from center
    # at angles 30° + k*60° for k=0..5
    # Hex circumradius = center-to-vertex distance
    # For hex with center-to-center spacing d = √3*s, circumradius = s (the size param)
    # But our axial size param s_level is different from the hex circumradius.
    # Hex circumradius = s_level (center-to-vertex = size parameter)
    # Actually for a regular hex with side length a:
    #   circumradius (center-to-vertex) = a
    #   center-to-center distance = √3 * a
    # Our s_level is the "size parameter" used in axial_to_cart.
    # In axial_to_cart, (1,0) maps to (√3*s, 0) → distance √3*s from origin
    # So center-to-center = √3*s, meaning the hex has side = s and circumradius = s.
    hex_R = s_res
    total_rotation = -res * IGEO7_THETA7

    coords = Vector{LonLat{Float64}}(undef, n_verts + 1)
    for k in 0:(n_verts - 1)
        # Pointy-top hex: first vertex at 30°
        angle = deg2rad(30.0 + k * (360.0 / n_verts)) + total_rotation
        vx = cx + hex_R * cos(angle)
        vy = cy + hex_R * sin(angle)
        p = _lambert_inverse(vx, vy, base1)
        coords[k + 1] = _xyz_to_lonlat(p)
    end
    coords[end] = coords[1]  # close polygon
    return coords
end

"""True if the cell is a pentagon (all digits are 0)."""
function _igeo7_is_pentagon(idx::UInt64)
    res = igeo7_resolution(idx)
    for i in 1:res
        igeo7_digit(idx, i) != 0 && return false
    end
    return true
end

#-----------------------------------------------------------------------------# AbstractCell Interface
grid(::IGEO7Cell) = IGEO7Grid()

resolution(o::IGEO7Cell) = igeo7_resolution(o.index)

is_pentagon(o::IGEO7Cell) = _igeo7_is_pentagon(o.index)

icon(o::IGEO7Cell) = is_pentagon(o) ? styled"{bright_red: ⬠}" : styled"{bright_green: ⬡}"

decode(o::IGEO7Cell) = string(igeo7_base_cell(o.index), "-", join(igeo7_digits(o.index)))
digits(o::IGEO7Cell) = igeo7_digits(o.index)
base_cell(o::IGEO7Cell) = igeo7_base_cell(o.index)
encode(o::IGEO7Cell) = igeo7_string(o.index)

#-----------------------------------------------------------------------------# GeoInterface
function GI.centroid(::GI.PolygonTrait, o::IGEO7Cell)
    return _igeo7_centroid(o.index)
end

function GI.coordinates(::GI.PolygonTrait, o::IGEO7Cell)
    return _igeo7_boundary(o.index)
end

function GI.area(::GI.PolygonTrait, o::IGEO7Cell)
    # Equal-area: total sphere surface / n_cells
    r = igeo7_resolution(o.index)
    n = 10 * 7.0^r + 2
    return 4π * EARTH_RADIUS_AUTHALIC^2 / n
end

#-----------------------------------------------------------------------------# Navigation
function parent(o::IGEO7Cell)
    r = resolution(o)
    r == 0 && return nothing
    return IGEO7Cell(igeo7_base_cell(o.index), igeo7_digits(o.index)[1:end-1])
end

function children(o::IGEO7Cell)
    r = resolution(o)
    r >= IGEO7_MAX_RES && return IGEO7Cell[]
    base = igeo7_base_cell(o.index)
    ds = igeo7_digits(o.index)
    pent = is_pentagon(o)
    n_children = pent ? 6 : 7
    return [IGEO7Cell(base, [ds; d]) for d in 0:(n_children - 1)]
end

function siblings(o::IGEO7Cell)
    r = resolution(o)
    r == 0 && return nothing
    return filter!(!=(o), children(parent(o)))
end

#-----------------------------------------------------------------------------# Misc
igeo7_pentagons(r::Integer) = [IGEO7Cell(b, fill(0, r)) for b in 0:11]
pentagons(::IGEO7Grid, r::Integer) = igeo7_pentagons(r)
ncells(::IGEO7Grid, res::Integer) = igeo7_n_cells(res)

#-----------------------------------------------------------------------------# igeo7cells
"""
    igeo7cells(geom, res; kw...)

Convenience wrapper for `cells(IGEO7Cell, geom, res; kw...)`.
"""
igeo7cells(args...; kw...) = cells(IGEO7Cell, args...; kw...)

function cells(::Type{IGEO7Cell}, ::GI.LineTrait, geom, res::Integer; containment=:center)
    containment in (:center, :shortest_path) ||
        throw(ArgumentError("Invalid containment mode for IGEO7 LineTrait. Expected :center or :shortest_path. Found: $containment."))
    coords = [LonLat(x...) for x in GI.coordinates(geom)]
    a_ll, b_ll = coords[1], coords[2]
    # Densify the line: sample points along the great circle and collect unique cells
    d = haversine(a_ll, b_ll)
    n_samples = max(2, ceil(Int, d / 1000.0))  # ~1 km spacing, at least 2 points
    az = azimuth(a_ll, b_ll)
    out = Set{IGEO7Cell}()
    for i in 0:n_samples
        frac = i / n_samples
        pt = destination(a_ll, az, frac * d)
        push!(out, IGEO7Cell(pt, res))
    end
    return collect(out)
end

function cells(::Type{IGEO7Cell}, ::GI.PolygonTrait, geom, res::Integer; containment=:center)
    containment in (:center, :overlap) ||
        throw(ArgumentError("Invalid containment mode for IGEO7 PolygonTrait. Expected :center or :overlap. Found: $containment."))
    ext = GI.extent(geom)
    X, Y = ext.X, ext.Y
    # Estimate step size from cell area: use 1/3 of cell diameter for dense sampling
    n = 10 * 7.0^res + 2
    cell_area_rad2 = 4π / n
    step_deg = rad2deg(sqrt(cell_area_rad2)) / 3.0
    step_deg = max(step_deg, 1e-6)
    out = Set{IGEO7Cell}()
    lat = Y[1]
    while lat <= Y[2]
        lon = X[1]
        while lon <= X[2]
            c = IGEO7Cell(LonLat(lon, lat), res)
            ct = GI.centroid(c)
            if GO.contains(geom, ct)
                push!(out, c)
            end
            lon += step_deg
        end
        lat += step_deg
    end
    return collect(out)
end
