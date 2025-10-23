module IGEO7

import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import GeometryBasics as GB
import ..GlobalGrids: Icosahedron
import StyledStrings: @styled_str

#-----------------------------------------------------------------------------# Grid
struct Grid{T}
    ico::Icosahedron{T}
    Grid{T}() where {T} = new{T}(Icosahedron{T}())
end
Grid() = Grid{Float64}()

Base.show(io::IO, o::Grid{T}) where {T} = print(io, "IGEO7.Grid{$T}")

# Center is (0, 0) in Axial Coordinates
struct Cell{T}
    grid::Grid{T}
    base::Int
    digits::Vector{Int}
    center::GB.Point3{T}
end
function Base.show(io::IO, o::Cell{T}) where {T}
    print(io, "Cell{$T}", styled"{bright_green: base=$(o.base)} {bright_cyan:digits=$(o.digits)}")
end


const scale_factor = 1 / sqrt(7)

function hexagon_axial_centers(o::Cell{T}) where {T}
    a = GB.Point2{T}(1, 0)
    b = GB.Point2{T}(0.5, sqrt(T(3)) / T(2))
    scale_factor .* (GB.Point2{T}(0, 0), a, (a-b), -b, -a, (b-a), b)
end

Base.getindex(grid::Grid{T}, i::Int) where {T} = Cell{T}(grid, i, Int[])
Base.getindex(cell::Cell{T}, digit::Int) where {T} = Cell{T}(cell.grid, cell.base, vcat(cell.digits, digit))


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


end
