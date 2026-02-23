#-----------------------------------------------------------------------------# Bit Layout
# Aperture 3/4: base(4 bits) + reserved(4 bits) + 28×2-bit digits
# Aperture 7/43: base(4 bits) + 20×3-bit digits

const DGG_BASE_OFFSET = 60  # bits [63:60] for base vertex (0-11)

# Per-aperture compile-time constants
@inline _dgg_bpd(::Val{3}) = 2;    @inline _dgg_bpd(::Val{4}) = 2
@inline _dgg_bpd(::Val{7}) = 3;    @inline _dgg_bpd(::Val{43}) = 3

@inline _dgg_maxd(::Val{3}) = 28;   @inline _dgg_maxd(::Val{4}) = 28
@inline _dgg_maxd(::Val{7}) = 20;   @inline _dgg_maxd(::Val{43}) = 20

@inline _dgg_mask(::Val{3}) = UInt64(0x3);  @inline _dgg_mask(::Val{4}) = UInt64(0x3)
@inline _dgg_mask(::Val{7}) = UInt64(0x7);  @inline _dgg_mask(::Val{43}) = UInt64(0x7)

@inline _dgg_maxval(::Val{3}) = 2;  @inline _dgg_maxval(::Val{4}) = 3
@inline _dgg_maxval(::Val{7}) = 6;  @inline _dgg_maxval(::Val{43}) = 6

#-----------------------------------------------------------------------------# DGGGrid
"""
    DGGGrid{Proj, Ap, Topo}

A discrete global grid parameterized by projection (`Proj`), aperture (`Ap`),
and topology (`Topo`).

- `Proj`: `:isea` (Snyder equal-area) or `:fuller` (Dymaxion)
- `Ap`: `3`, `4`, `7`, or `43` (mixed 4-3)
- `Topo`: `:hex`, `:tri`, or `:diamond`

### Examples

    ISEA3H()   # ISEA aperture-3 hexagonal grid
"""
struct DGGGrid{Proj, Ap, Topo} <: AbstractGrid end

GI.crs(::DGGGrid) = GFT.EPSG(4326)
GI.ncoord(::DGGGrid) = 2

#-----------------------------------------------------------------------------# DGGCell
"""
    DGGCell{Proj, Ap, Topo}(index::UInt64)
    DGGCell{Proj, Ap, Topo}(base::Integer, digits::AbstractVector{<:Integer})
    DGGCell{Proj, Ap, Topo}(coord::LonLat, res::Integer)
    DGGCell{Proj, Ap, Topo}(str::AbstractString)

A single cell in a DGG grid.  Wraps a `UInt64` index encoding a 4-bit base
cell (0–11) and hierarchical refinement digits.

### Examples

    ISEA3HCell(LonLat(-75.0, 54.0), 5)
    ISEA3HCell(0, [1, 2, 0])
"""
struct DGGCell{Proj, Ap, Topo} <: AbstractCell
    index::UInt64
end

#-----------------------------------------------------------------------------# Grid/Cell Aliases
const ISEA3H    = DGGGrid{:isea, 3, :hex}
const ISEA4H    = DGGGrid{:isea, 4, :hex}
const ISEA7H    = DGGGrid{:isea, 7, :hex}
const ISEA43H   = DGGGrid{:isea, 43, :hex}
const ISEA4T    = DGGGrid{:isea, 4, :tri}
const ISEA4D    = DGGGrid{:isea, 4, :diamond}
const FULLER3H  = DGGGrid{:fuller, 3, :hex}
const FULLER4H  = DGGGrid{:fuller, 4, :hex}
const FULLER7H  = DGGGrid{:fuller, 7, :hex}
const FULLER43H = DGGGrid{:fuller, 43, :hex}
const FULLER4T  = DGGGrid{:fuller, 4, :tri}
const FULLER4D  = DGGGrid{:fuller, 4, :diamond}

const ISEA3HCell    = DGGCell{:isea, 3, :hex}
const ISEA4HCell    = DGGCell{:isea, 4, :hex}
const ISEA7HCell    = DGGCell{:isea, 7, :hex}
const ISEA43HCell   = DGGCell{:isea, 43, :hex}
const ISEA4TCell    = DGGCell{:isea, 4, :tri}
const ISEA4DCell    = DGGCell{:isea, 4, :diamond}
const FULLER3HCell  = DGGCell{:fuller, 3, :hex}
const FULLER4HCell  = DGGCell{:fuller, 4, :hex}
const FULLER7HCell  = DGGCell{:fuller, 7, :hex}
const FULLER43HCell = DGGCell{:fuller, 43, :hex}
const FULLER4TCell  = DGGCell{:fuller, 4, :tri}
const FULLER4DCell  = DGGCell{:fuller, 4, :diamond}

#-----------------------------------------------------------------------------# Bit Manipulation
"""Extract base vertex (0–11) from raw index."""
@inline dgg_base(idx::UInt64) = Int((idx >> DGG_BASE_OFFSET) & 0xF)

"""Get digit at 1-based position `i`."""
@inline function dgg_digit(idx::UInt64, i::Int, ::Val{Ap}) where Ap
    bpd = _dgg_bpd(Val(Ap))
    maxd = _dgg_maxd(Val(Ap))
    shift = bpd * (maxd - i)
    Int((idx >> shift) & _dgg_mask(Val(Ap)))
end

"""Resolution = number of leading non-padding digits."""
function dgg_resolution(idx::UInt64, ::Val{Ap}) where Ap
    pad = Int(_dgg_mask(Val(Ap)))
    maxd = _dgg_maxd(Val(Ap))
    for i in 1:maxd
        dgg_digit(idx, i, Val(Ap)) == pad && return i - 1
    end
    return maxd
end

"""Return active digits as a vector."""
function dgg_digits(idx::UInt64, ::Val{Ap}) where Ap
    r = dgg_resolution(idx, Val(Ap))
    [dgg_digit(idx, i, Val(Ap)) for i in 1:r]
end

"""Hex string of raw index."""
dgg_string(idx::UInt64) = string(idx, base=16)
dgg_string(o::DGGCell) = dgg_string(o.index)

"""Number of cells at resolution `r` for hex grid with aperture `a`."""
dgg_n_cells(a::Integer, r::Integer) = 10 * big(a)^r + 2

#-----------------------------------------------------------------------------# Constructors
function DGGCell{P,A,T}(base::Integer, digits::AbstractVector{<:Integer}) where {P,A,T}
    0 <= base <= 11 || throw(ArgumentError("base must be in 0:11, got $base"))
    maxd = _dgg_maxd(Val(A))
    mask = _dgg_mask(Val(A))
    maxval = _dgg_maxval(Val(A))
    bpd = _dgg_bpd(Val(A))
    length(digits) <= maxd || throw(ArgumentError("digits length must be <= $maxd"))
    all(d -> 0 <= d <= maxval, digits) || throw(ArgumentError("each digit must be in 0:$maxval"))
    idx = UInt64(base) << DGG_BASE_OFFSET
    for i in 1:maxd
        d = i <= length(digits) ? UInt64(digits[i]) : mask  # pad with all-1s
        shift = bpd * (maxd - i)
        idx |= (d & mask) << shift
    end
    return DGGCell{P,A,T}(idx)
end

function DGGCell{P,A,T}(str::AbstractString) where {P,A,T}
    DGGCell{P,A,T}(parse(UInt64, str; base=16))
end

#-----------------------------------------------------------------------------# AbstractCell Interface
grid(::DGGCell{P,A,T}) where {P,A,T} = DGGGrid{P,A,T}()

resolution(o::DGGCell{P,A,T}) where {P,A,T} = dgg_resolution(o.index, Val(A))

function _dgg_is_pentagon(idx::UInt64, ::Val{Ap}) where Ap
    r = dgg_resolution(idx, Val(Ap))
    for i in 1:r
        dgg_digit(idx, i, Val(Ap)) != 0 && return false
    end
    return true
end

is_pentagon(o::DGGCell{P,A,T}) where {P,A,T} = _dgg_is_pentagon(o.index, Val(A))

icon(o::DGGCell) = is_pentagon(o) ? styled"{bright_red: ⬠}" : styled"{bright_green: ⬡}"

decode(o::DGGCell{P,A,T}) where {P,A,T} = string(dgg_base(o.index), "-", join(dgg_digits(o.index, Val(A))))
digits(o::DGGCell{P,A,T}) where {P,A,T} = dgg_digits(o.index, Val(A))
base_cell(o::DGGCell) = dgg_base(o.index)
ncells(::DGGGrid{P,A,:hex}, res::Integer) where {P,A} = dgg_n_cells(A, res)
encode(o::DGGCell) = dgg_string(o.index)
