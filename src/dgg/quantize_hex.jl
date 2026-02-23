#-----------------------------------------------------------------------------# Aperture-3 Hex Constants

# Rotation angle per level: 30° (Class I ↔ Class III alternation)
const DGG_AP3_THETA = π / 6

# Digit offsets in axial (q, r) coordinates — coset representatives of Z²/M₃Z²
const DGG_AP3_DIGIT_OFFSETS = SVector{3, SVector{2, Int}}(
    SVector(0, 0),   # digit 0: center
    SVector(1, 0),   # digit 1
    SVector(0, 1),   # digit 2
)

# Hex spacing at resolution 1 for aperture 3.
# Constrained so that level-1 child centers (at distance √3*S1 from base vertex)
# stay within the vertex Voronoi region (apothem ≈ 0.547 in Lambert coords).
const DGG_AP3_S1 = let
    v1 = IGEO7_VERTICES_3D[1]
    v2 = IGEO7_VERTICES_3D[IGEO7_VERTEX_NEIGHBORS[1][1]]
    edge_angle = acos(clamp(dot(v1, v2), -1.0, 1.0))
    apothem = 2.0 * sin(edge_angle / 4.0)
    apothem / sqrt(3.0) * 0.95
end

#-----------------------------------------------------------------------------# Aperture-specific accessors (extend for aperture 4, 7 later)
_dgg_theta(::Val{3})          = DGG_AP3_THETA
_dgg_digit_offsets(::Val{3})  = DGG_AP3_DIGIT_OFFSETS
_dgg_s1(::Val{3})             = DGG_AP3_S1
_dgg_scale_factor(::Val{3})   = sqrt(3.0)
_dgg_n_children(::Val{3})     = 3

#-----------------------------------------------------------------------------# Bottom-up digit extraction (aperture-3 matrix algebra)

# M₃ = [2 1; -1 1] — maps coarser-level axial to finer-level axial (det=3, scale √3, rotate -π/6)
@inline _dgg_apply_M(q::Int, r::Int, ::Val{3}) = (2q + r, -q + r)

# M₃⁻¹ = [1 -1; 1 2] / 3
@inline _dgg_apply_M_inv(q::Int, r::Int, ::Val{3}) = (div(q - r, 3), div(q + 2r, 3))

# Digit = coset index: mod(q - r, 3) for aperture 3
@inline _dgg_axial_to_digit(q::Int, r::Int, ::Val{3}) = mod(q - r, 3)

"""Extract hierarchical digits bottom-up from finest-level axial coordinates."""
function _dgg_bottomup_digits(q::Int, r::Int, res::Integer, ::Val{A}) where A
    digits = Vector{Int}(undef, res)
    for level in res:-1:1
        d = _dgg_axial_to_digit(q, r, Val(A))
        digits[level] = d
        off = _dgg_digit_offsets(Val(A))[d + 1]
        q, r = _dgg_apply_M_inv(q - off[1], r - off[2], Val(A))
    end
    return digits
end

"""Reconstruct finest-level axial coordinates from cell index (top-down digit application)."""
function _dgg_finest_hex(idx::UInt64, ::Val{A}) where A
    base1 = dgg_base(idx) + 1
    res = dgg_resolution(idx, Val(A))
    q, r = 0, 0
    for level in 1:res
        d = dgg_digit(idx, level, Val(A))
        off = _dgg_digit_offsets(Val(A))[d + 1]
        nq, nr = _dgg_apply_M(q, r, Val(A))
        q = nq + off[1]
        r = nr + off[2]
    end
    return base1, q, r
end
