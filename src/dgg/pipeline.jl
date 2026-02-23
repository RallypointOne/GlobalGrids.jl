#-----------------------------------------------------------------------------# Forward Transform: LonLat -> DGGCell (hex grids)

function DGGCell{P,A,:hex}(coord::LonLat, res::Integer) where {P,A}
    maxd = _dgg_maxd(Val(A))
    0 <= res <= maxd || throw(ArgumentError("resolution must be in 0:$maxd, got $res"))
    _dgg_forward_hex(coord, res, Val(P), Val(A))
end

DGGCell{P,A,:hex}((x, y)::NTuple{2, Real}, res::Integer) where {P,A} = DGGCell{P,A,:hex}(LonLat(Float64(x), Float64(y)), res)

"""Compute a candidate DGG cell on face `fi` at resolution `res` from Lambert coords.
Returns the cell index (UInt64)."""
@inline function _dgg_quantize_on_face(px::Float64, py::Float64, fi::Int, res::Integer, ::Val{A}) where A
    s1 = _dgg_s1(Val(A))
    scale = _dgg_scale_factor(Val(A))
    theta = _dgg_theta(Val(A))

    s_finest = s1 / scale^(res - 1)
    total_rotation = res * theta

    # Undo cumulative rotation to align with finest-level hex grid
    rx, ry = _rotate2d(px, py, -total_rotation)

    # Convert to fractional axial, snap to nearest hex center
    fq, fr = _cart_to_axial(rx, ry, s_finest)
    q, r = _hex_round(fq, fr)

    # Extract digits bottom-up using modular arithmetic
    digits = _dgg_bottomup_digits(q, r, res, Val(A))

    return DGGCell{:isea,A,:hex}(fi - 1, digits).index
end

"""Forward transform: quantize a LonLat to the nearest DGG hex cell.

Uses face-centered Lambert projection with bottom-up digit extraction.
Tries multiple candidate faces near the input point and picks the cell
whose centroid is closest on the sphere, handling face boundary cases."""
function _dgg_forward_hex(coord::LonLat, res::Integer, ::Val{P}, ::Val{A}) where {P,A}
    p = _lonlat_to_xyz(coord)

    if res == 0
        fi = _nearest_face(p)
        return DGGCell{P,A,:hex}(fi - 1, Int[])
    end

    # Score all faces by dot product with input point, try top candidates
    best_idx = UInt64(0)
    best_dot = -Inf

    for fi in 1:20
        face_dot = dot(p, DGG_FACE_CENTROIDS[fi])
        # Skip faces whose centroid is too far (cos > 0.3 means < ~72 degrees)
        face_dot < 0.3 && continue

        px, py = _face_lambert_forward(p, fi)
        idx = _dgg_quantize_on_face(px, py, fi, res, Val(A))

        # Score: dot product of cell centroid with input point (higher = closer)
        centroid_p = _lonlat_to_xyz(_dgg_centroid(idx, Val(A)))
        cdot = dot(p, centroid_p)
        if cdot > best_dot
            best_dot = cdot
            best_idx = idx
        end
    end

    return DGGCell{P,A,:hex}(best_idx)
end

#-----------------------------------------------------------------------------# Inverse Transform: DGGCell -> center LonLat

"""Compute face Lambert-plane center (x, y) and face index (1-indexed) for a hex DGG cell.

Accumulates rotated offsets top-down, matching the forward transform's geometric model."""
function _dgg_center_lambert(idx::UInt64, ::Val{A}) where A
    fi = dgg_base(idx) + 1  # face index, 1-indexed
    res = dgg_resolution(idx, Val(A))
    res == 0 && return fi, 0.0, 0.0
    offsets = _dgg_digit_offsets(Val(A))
    s1 = _dgg_s1(Val(A))
    scale = _dgg_scale_factor(Val(A))
    theta = _dgg_theta(Val(A))
    cx, cy = 0.0, 0.0
    for level in 1:res
        d = dgg_digit(idx, level, Val(A))
        s_level = s1 / scale^(level - 1)
        rotation = level * theta
        off = offsets[d + 1]
        ox, oy = _axial_to_cart(off[1], off[2], s_level)
        ox, oy = _rotate2d(ox, oy, rotation)
        cx += ox
        cy += oy
    end
    return fi, cx, cy
end

function _dgg_centroid(idx::UInt64, ::Val{A}) where A
    fi, cx, cy = _dgg_center_lambert(idx, Val(A))
    p = _face_lambert_inverse(cx, cy, fi)
    return _xyz_to_lonlat(p)
end
