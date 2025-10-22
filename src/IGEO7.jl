module IGEO7

import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import GeometryBasics as GB
import ..GlobalGrids: Icosahedron, radius, R
import StyledStrings: @styled_str

#-----------------------------------------------------------------------------# Grid
struct Grid{T}
    ico::Icosahedron{T}
    function Grid{T}(radius=R) where {T}
        ico = Icosahedron{T}(radius)
        return new{T}(ico)
    end
end
Grid(radius=R) = Grid{Float64}(radius)
radius(Grid::Grid) = radius(Grid.ico)

function Base.show(io::IO, o::Grid)
    print(io, styled"{bright_green:IGEO7.Grid:} $(o.ico)")
end


#-----------------------------------------------------------------------------# Z7 Indexing
# Constants from IGEO7 Z7 spec
const Z7_BASE_BITS  = 4
const Z7_DIGIT_BITS = 3
const Z7_MAX_RES    = 20
const Z7_PAD_VALUE  = 0x07

"""
    decode(idx::UInt64)

Split a Z7 index into its parts:
- `base`: base cell id (0–11)
- `digits`: Vector{UInt8} of up to 20 base-7 digits (values 0–6)
- `resolution`: number of valid digits before padding
- `padding`: number of padded digits (= 7)

Assumes 64-bit Z7 format:
[ 4 bits base ][ 20 * 3-bit digits (LSB) ]
"""
function decode(idx::UInt64)
    # Extract base cell (top 4 bits)
    base = idx >> (Z7_MAX_RES * Z7_DIGIT_BITS)

    # Extract 20 refinement digits (3 bits each)
    digits = UInt8[]
    mask = (UInt64(1) << Z7_DIGIT_BITS) - 1  # 0b111
    for i in 0:(Z7_MAX_RES-1)
        d = UInt8((idx >> (i * Z7_DIGIT_BITS)) & mask)
        push!(digits, d)
    end

    # Resolution = number of digits before first padding value (7)
    resolution = findfirst(==(Z7_PAD_VALUE), digits)
    resolution = isnothing(resolution) ? Z7_MAX_RES : resolution - 1
    padding = Z7_MAX_RES - resolution

    return (base = UInt8(base), digits = digits, resolution = resolution, padding = padding)
end

"""
    encode(base_cell, digits)

Construct a Z7 index from:
- `base_cell`: Integer in 0-11.
- `digits`: collection (maximum length 20) of Integers 0-6.
"""
function encode(base_cell, digits)
    -1 < base_cell < 12 || throw(ArgumentError("base_cell must be in 0-11.  Found: $base_cell"))
    all(x -> (-1 < x < 7), digits) || throw(ArgumentError("digits must all be in 0-6.  Found: $digits"))
    idx = UInt64(base_cell - 1) << (Z7_MAX_RES * Z7_DIGIT_BITS)
    for (i, d) in enumerate(digits)
        idx |= (UInt64(d) << ((i - 1) * Z7_DIGIT_BITS))
    end
    for i in (length(digits)+1):Z7_MAX_RES
        idx |= (UInt64(Z7_PAD_VALUE) << ((i - 1) * Z7_DIGIT_BITS))
    end
    return idx
end

function z7string(x::UInt64)
    (; base, digits) = decode(x)
    filter!(!=(0x07), digits)
    string(lpad(Int(base), 2, '0'), join(Int.(digits), ""))
end

#-----------------------------------------------------------------------------# Cell
struct Cell
    index::UInt64  # Z7 index
end



end
