#-----------------------------------------------------------------------------# Cell Boundary (polygon vertices) for hex DGG cells

function _dgg_boundary(idx::UInt64, ::Val{A}) where A
    fi, cx, cy = _dgg_center_lambert(idx, Val(A))
    res = dgg_resolution(idx, Val(A))
    n_verts = 6

    s1 = _dgg_s1(Val(A))
    scale = _dgg_scale_factor(Val(A))
    theta = _dgg_theta(Val(A))
    s_res = res == 0 ? s1 : s1 / scale^(res - 1)
    hex_R = s_res
    total_rotation = res * theta

    coords = Vector{LonLat{Float64}}(undef, n_verts + 1)
    for k in 0:(n_verts - 1)
        angle = deg2rad(30.0 + k * (360.0 / n_verts)) + total_rotation
        vx = cx + hex_R * cos(angle)
        vy = cy + hex_R * sin(angle)
        p = _face_lambert_inverse(vx, vy, fi)
        coords[k + 1] = _xyz_to_lonlat(p)
    end
    coords[end] = coords[1]  # close polygon
    return coords
end

#-----------------------------------------------------------------------------# GeoInterface
GI.centroid(::GI.PolygonTrait, o::DGGCell{P,A,T}) where {P,A,T} = _dgg_centroid(o.index, Val(A))

GI.coordinates(::GI.PolygonTrait, o::DGGCell{P,A,T}) where {P,A,T} = _dgg_boundary(o.index, Val(A))

function GI.area(::GI.PolygonTrait, o::DGGCell{P,A,:hex}) where {P,A}
    r = resolution(o)
    n = 20 * Float64(A)^r
    return 4π * EARTH_RADIUS_AUTHALIC^2 / n
end
