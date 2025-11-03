module LibH3

import ..GlobalGrids: AbstractCoordinates, LonLatAuthalic
import GeoInterface as GI
import H3_jll: libh3
import GeometryBasics: Point

#-----------------------------------------------------------------------------# API
const H3Error = UInt32
const H3Index = UInt64

#-----------------------------------------------------------------------------# Types
# Internal type for H3 latitude/longitude coordinates in *radians*.  Not intended for public use.
struct LatLng
    lat::Cdouble
    lng::Cdouble
end
LatLng(o::AbstractCoordinates) = LatLng(LonLatAuthalic(o))
LatLng(o::LonLatAuthalic) = LatLng(deg2rad(o[2]), deg2rad(o[1]))  # (lon, lat) in degrees
LonLatAuthalic(o::LatLng) = LonLatAuthalic(Point(rad2deg(o.lng), rad2deg(o.lat)))  # (lon, lat) in degrees

# H3 Cell vertices
struct CellBoundary
    numVerts::Cint
    verts::NTuple{10, LatLng}
end

# Note: verts[1] does not need to be repeated at verts[end]
struct GeoLoop
    numVerts::Cint
    verts::Ptr{LatLng}
end

struct GeoPolygon
    geoloop::GeoLoop
    numHoles::Cint
    holes::Ptr{GeoLoop}
end

#-----------------------------------------------------------------------------# Functions
macro libh3(ex)
    fun, ret = ex.args
    fun_name, args... = fun.args
    arg_names = [x.args[1] for x in args]
    arg_types = Expr(:tuple, [x.args[2] for x in args]...)
    ret_check_ex = ret == :H3Error ? :(if res != 0
            throw(ErrorException("H3 call to $(string($fun_name)) failed with error code $res."))
        end) : nothing
    esc(quote
        function $fun_name($(arg_names...))
            res = ccall(($(QuoteNode(fun_name)), libh3), $ret, $arg_types, $(arg_names...))
            $ret_check_ex
            return res
        end
    end)
end

@libh3 latLngToCell(g::Ptr{LatLng}, res::Cint, out::Ptr{UInt64})::H3Error
@libh3 cellToLatLng(cell::H3Index, g::Ptr{LatLng})::H3Error
@libh3 cellToBoundary(cell::H3Index, bndry::Ptr{CellBoundary})::H3Error
@libh3 isValidCell(cell::H3Index)::Cint
@libh3 isPentagon(cell::H3Index)::Cint
@libh3 getBaseCellNumber(cell::H3Index)::Cint
@libh3 getResolution(cell::H3Index)::Cint
@libh3 stringToH3(str::Cstring, out::Ptr{H3Index})::H3Error
@libh3 maxFaceCount(h3::H3Index, out::Ptr{Cint})::H3Error
@libh3 getIcosahedronFaces(h::H3Index, out::Ptr{Cint})::H3Error
@libh3 gridDistance(origin::H3Index, h3::H3Index, distance::Ptr{Int64})::H3Error
@libh3 gridRing(origin::H3Index, k::Cint, out::Ptr{H3Index})::H3Error
@libh3 gridDisk(origin::H3Index, k::Cint, out::Ptr{H3Index})::H3Error
@libh3 gridPathCells(start::H3Index, _end::H3Index, out::Ptr{H3Index})::H3Error
@libh3 gridPathCellsSize(start::H3Index, _end::H3Index, size::Ptr{Int64})::H3Error
@libh3 cellToChildren(cell::H3Index, childRes::Cint, children::Ptr{H3Index})::H3Error
@libh3 cellToChildrenSize(cell::H3Index, childRes::Cint, out::Ptr{Int64})::H3Error

@libh3 polygonToCellsExperimental(geoPolygon::Ptr{GeoPolygon}, res::Cint, flags::UInt32, size::Int64, out::Ptr{H3Index})::H3Error
@libh3 maxPolygonToCellsSizeExperimental(geoPolygon::Ptr{GeoPolygon}, res::Cint, flags::UInt32, out::Ptr{Int64})::H3Error


end  # module LibH3
