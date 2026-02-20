#-----------------------------------------------------------------------------# rotation_matrix
"""
    rotation_matrix(axis::SVector{3}, angle)

Create a 3x3 rotation matrix using Rodrigues' rotation formula.
Rotates by `angle` (radians) around `axis` (must be normalized).
"""
function rotation_matrix(axis::SVector{3,T}, angle) where {T}
    K = SMatrix{3,3,T}(0, axis[3], -axis[2], -axis[3], 0, axis[1], axis[2], -axis[1], 0)
    return I + sin(angle) * K + (1 - cos(angle)) * (K * K)
end

#-----------------------------------------------------------------------------# icosahedron
# Base icosahedron with circumradius = 1
function unit_icosahedron(T = Float64)
    φ = (1 + sqrt(5)) / 2
    verts = normalize.(Point{3,T}[
        (-1,  φ, 0),  (1,  φ, 0),  (-1, -φ, 0),  (1, -φ, 0), (0, -1,  φ), (0,  1,  φ),
        (0, -1, -φ), (0,  1, -φ), (φ,  0, -1), (φ,  0,  1), (-φ, 0, -1), (-φ, 0,  1)
    ])
    faces = GB.TriangleFace.([(1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12), (2,6,10), (6,12,5),
        (12,11,3), (11,8,7), (8,2,9), (4,10,5), (4,5,3), (4,3,7), (4,7,9),  (4,9,10), (5,10,6),
        (3,5,12), (7,3,11), (9,7,8), (10,9,2)
    ])
    return Mesh(verts, faces)
end

"""
    icosahedron(vertex=Point3d(0,0,1), azimuth=0.0)
    icosahedron(preset::Symbol)

Create an icosahedron mesh on the unit sphere.  A vertex is placed at `vertex`
and the mesh is rotated by `azimuth` (radians) around that axis.

Preset orientations: `:poles`, `:isea`, `:dymaxion`.

### Examples

    icosahedron()
    icosahedron(:isea)
"""
function icosahedron(vertex::Point3{T} = Point3d(0,0,1), azimuth::T = 0.0) where {T}
    unit_ico = unit_icosahedron(eltype(vertex))
    ref_vertex = normalize(SVector{3,T}(unit_ico.position[1]...))  # First vertex of unit_ico
    target_dir = normalize(SVector{3,T}(vertex...))  # Target direction and magnitude

    # Rotation 1: Align ref_vertex with target direction
    axis1 = cross(ref_vertex, target_dir)
    if norm(axis1) < eps(T) * 10
        # Vertices are parallel or antiparallel
        if dot(ref_vertex, target_dir) > 0
            R1 = SMatrix{3,3,T}(I)  # No rotation needed
        else
            # 180° rotation - pick perpendicular axis
            perp = abs(ref_vertex[1]) < T(0.9) ? SVector{3,T}(1, 0, 0) : SVector{3,T}(0, 1, 0)
            axis1 = normalize(cross(ref_vertex, perp))
            R1 = rotation_matrix(axis1, T(π))
        end
    else
        axis1 = normalize(axis1)
        angle1 = acos(clamp(dot(ref_vertex, target_dir), T(-1), T(1)))
        R1 = rotation_matrix(axis1, angle1)
    end

    # Rotation 2: Rotate around the vertex axis by azimuth
    R2 = rotation_matrix(target_dir, azimuth)

    Rot = R2 * R1  # Combined transformation: scale and rotate

    new_vertices = normalize.([Point{3,T}((Rot * SVector{3,T}(v...))...) for v in unit_ico.position])

    return Mesh(new_vertices, unit_ico.faces)
end

function icosahedron(preset::Symbol)
    preset == :poles && return icosahedron(ORIENTATION.POLES...)
    preset == :isea  && return icosahedron(ORIENTATION.ISEA...)
    preset == :dymaxion && return icosahedron(ORIENTATION.DYMAXION...)
    throw(ArgumentError("Unknown preset orientation: $preset"))
end

unit_sphere(T = Float64) = GB.Sphere(Point3{T}(0.0, 0.0, 0.0), 0.0)

#-----------------------------------------------------------------------------# preset orientations
# Reference: https://discreteglobalgrids.org/dgg-orientation/
const ORIENTATION = (
    POLES = (Point3d(0,0, 1), 0.0),
    ISEA = (Point3d(0.518136909297074, 0.1030638392550601, 0.8490653615959627), 0.0),
    DYMAXION = (Point3d(0.9950201403480873, -0.09134877248651022, 0.039878842346294276), deg2rad(7.466580))  # LonLat(-5.245390, 2.3008820)
)


#-----------------------------------------------------------------------------# geometry
# Dual of icosahedron mesh: pentagonal faces at vertices
function dual(m::Mesh)
    verts = Point.(values(GB.face_normals(m.position, m.faces)))
    face_normals = m.position
    faces = GB.NgonFace{5,Int}[]
    map(face_normals) do n
        face = sort(verts, by = v -> -dot(v, n))[1:5]
        counterclockwise!(face)
        idx = [findfirst(==(v), verts) for v in face]
        push!(faces, GB.NgonFace(idx...))
    end
    return Mesh(verts, faces)
end

# Force collection of points to be in counter-clockwise order
function counterclockwise!(points::AbstractVector{<:Point})
    c = sum(points) / length(points)
    n = normalize(c)

    # Build local 2D coordinate system (u, v) in polygon plane
    tmp = abs(n[1]) < 0.9 ? SVector(1.0, 0.0, 0.0) : SVector(0.0, 1.0, 0.0)
    u = normalize(cross(n, tmp))
    v = cross(n, u)

    # Sort by polar angle in the local 2D system
    angles = [atan(dot(pt - c, v), dot(pt - c, u)) for pt in points]
    idx = sortperm(angles)
    points .= points[idx]
end
