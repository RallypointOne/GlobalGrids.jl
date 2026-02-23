#-----------------------------------------------------------------------------# Forward Transform: LonLat → DGGCell (hex grids)

function DGGCell{P,A,:hex}(coord::LonLat, res::Integer) where {P,A}
    maxd = _dgg_maxd(Val(A))
    0 <= res <= maxd || throw(ArgumentError("resolution must be in 0:$maxd, got $res"))
    _dgg_forward_hex(coord, res, Val(P), Val(A))
end

DGGCell{P,A,:hex}((x, y)::NTuple{2, Real}, res::Integer) where {P,A} = DGGCell{P,A,:hex}(LonLat(Float64(x), Float64(y)), res)

"""Find the 1-indexed base vertex nearest to 3D unit vector `p`."""
function _nearest_vertex(p::SVector{3,Float64})
    best_dot = -Inf
    base1 = 1
    for i in 1:12
        d = dot(p, IGEO7_VERTICES_3D[i])
        if d > best_dot
            best_dot = d
            base1 = i
        end
    end
    return base1
end

"""Quantize Lambert-plane point to finest-level hex and extract digits bottom-up."""
function _quantize_and_extract(px::Float64, py::Float64, base1::Int, res::Integer, ::Val{P}, ::Val{A}) where {P,A}
    s_finest = _dgg_s1(Val(A)) / _dgg_scale_factor(Val(A))^(res - 1)
    total_rot = -res * _dgg_theta(Val(A))

    # Undo total rotation to align with finest-level hex frame
    ux, uy = _rotate2d(px, py, -total_rot)

    # Convert to fractional axial and round to nearest hex
    fq, fr = _cart_to_axial(ux, uy, s_finest)
    q, r = _hex_round(fq, fr)

    # Bottom-up digit extraction
    digits = _dgg_bottomup_digits(q, r, res, Val(A))
    return DGGCell{P,A,:hex}(base1 - 1, digits)
end

function _dgg_forward_hex(coord::LonLat, res::Integer, ::Val{P}, ::Val{A}) where {P,A}
    p = _lonlat_to_xyz(coord)
    base1 = _nearest_vertex(p)

    res == 0 && return DGGCell{P,A,:hex}(base1 - 1, Int[])

    # Try nearest vertex first
    px, py = _lambert_forward(p, base1)
    best_cell = _quantize_and_extract(px, py, base1, res, Val(P), Val(A))
    best_centroid = _lonlat_to_xyz(_dgg_centroid(best_cell.index, Val(A)))
    best_dot = dot(p, best_centroid)

    # Search all other vertices for a cell whose centroid is closer to input.
    # This guarantees round-trips: forward(inverse(forward(x))) == forward(x)
    # because inverse(cell) produces centroid with dot=1.0 against itself.
    for v in 1:12
        v == base1 && continue
        px_v, py_v = _lambert_forward(p, v)
        cell_v = _quantize_and_extract(px_v, py_v, v, res, Val(P), Val(A))
        centroid_v = _lonlat_to_xyz(_dgg_centroid(cell_v.index, Val(A)))
        d = dot(p, centroid_v)
        if d > best_dot
            best_dot = d
            best_cell = cell_v
        end
    end

    return best_cell
end

#-----------------------------------------------------------------------------# Inverse Transform: DGGCell → center LonLat

"""Compute Lambert-plane center (x, y) and base vertex (1-indexed) for a hex DGG cell.

Uses bottom-up approach: convert digits to finest-level hex coords, then to Lambert.
This is consistent with the bottom-up forward transform, ensuring round-trip accuracy."""
function _dgg_center_lambert(idx::UInt64, ::Val{A}) where A
    base1, q, r = _dgg_finest_hex(idx, Val(A))
    res = dgg_resolution(idx, Val(A))
    res == 0 && return base1, 0.0, 0.0
    s_finest = _dgg_s1(Val(A)) / _dgg_scale_factor(Val(A))^(res - 1)
    total_rot = -res * _dgg_theta(Val(A))
    ux, uy = _axial_to_cart(q, r, s_finest)
    cx, cy = _rotate2d(ux, uy, total_rot)
    return base1, cx, cy
end

function _dgg_centroid(idx::UInt64, ::Val{A}) where A
    base1, cx, cy = _dgg_center_lambert(idx, Val(A))
    p = _lambert_inverse(cx, cy, base1)
    return _xyz_to_lonlat(p)
end

#-----------------------------------------------------------------------------# Geodesic helpers for DGGCell
haversine(a::DGGCell{P,A,T}, b::DGGCell{P,A,T}) where {P,A,T} = haversine(GI.centroid(a), GI.centroid(b))

function destination(a::DGGCell{P,A,T}, azimuth_deg, m) where {P,A,T}
    DGGCell{P,A,T}(destination(GI.centroid(a), azimuth_deg, m), resolution(a))
end
