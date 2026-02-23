#-----------------------------------------------------------------------------# dggcells convenience wrapper
"""
    dggcells(geom, res; grid=ISEA3H(), kw...)

Convenience wrapper for `cells(DGGCell{...}, geom, res; kw...)`.
"""
function dggcells(args...; grid::DGGGrid{P,A,T}=ISEA3H(), kw...) where {P,A,T}
    cells(DGGCell{P,A,T}, args...; kw...)
end

#-----------------------------------------------------------------------------# Line
function cells(::Type{DGGCell{P,A,:hex}}, ::GI.LineTrait, geom, res::Integer; containment=:center) where {P,A}
    containment in (:center, :shortest_path) ||
        throw(ArgumentError("Invalid containment mode. Expected :center or :shortest_path. Found: $containment."))
    coords = [LonLat(x...) for x in GI.coordinates(geom)]
    a_ll, b_ll = coords[1], coords[2]
    d = haversine(a_ll, b_ll)
    n_samples = max(2, ceil(Int, d / 1000.0))
    az = azimuth(a_ll, b_ll)
    out = Set{DGGCell{P,A,:hex}}()
    for i in 0:n_samples
        frac = i / n_samples
        pt = destination(a_ll, az, frac * d)
        push!(out, DGGCell{P,A,:hex}(pt, res))
    end
    return collect(out)
end

#-----------------------------------------------------------------------------# Polygon
function cells(::Type{DGGCell{P,A,:hex}}, ::GI.PolygonTrait, geom, res::Integer; containment=:center) where {P,A}
    containment in (:center, :overlap) ||
        throw(ArgumentError("Invalid containment mode. Expected :center or :overlap. Found: $containment."))
    ext = GI.extent(geom)
    X, Y = ext.X, ext.Y
    n = 10 * Float64(A)^res + 2
    cell_area_rad2 = 4π / n
    step_deg = rad2deg(sqrt(cell_area_rad2)) / 3.0
    step_deg = max(step_deg, 1e-6)
    out = Set{DGGCell{P,A,:hex}}()
    lat = Y[1]
    while lat <= Y[2]
        lon = X[1]
        while lon <= X[2]
            c = DGGCell{P,A,:hex}(LonLat(lon, lat), res)
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
