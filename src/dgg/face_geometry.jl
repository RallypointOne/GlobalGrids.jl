#-----------------------------------------------------------------------------# Icosahedron Face Geometry
# Precompute face centroids, vertex indices, and Lambert projection bases
# for face-centered DGG addressing.

"""Precompute ISEA icosahedron face centroids and vertex indices."""
function _compute_face_data()
    ico = icosahedron(:isea)
    verts_3d = [SVector{3,Float64}(v...) for v in ico.position]

    n_faces = length(ico.faces)
    centroids = Vector{SVector{3,Float64}}(undef, n_faces)
    vert_indices = Vector{SVector{3,Int}}(undef, n_faces)

    for (fi, face) in enumerate(ico.faces)
        i, j, k = face[1], face[2], face[3]
        vert_indices[fi] = SVector(i, j, k)
        centroids[fi] = normalize(verts_3d[i] + verts_3d[j] + verts_3d[k])
    end

    return SVector{20}(centroids), SVector{20}(vert_indices)
end

const DGG_FACE_CENTROIDS, DGG_FACE_VERTEX_INDICES = _compute_face_data()

"""Find icosahedron face (1-indexed) whose centroid is nearest to point `p`."""
function _nearest_face(p::SVector{3,Float64})
    best = -Inf
    best_fi = 1
    for fi in 1:20
        d = dot(p, DGG_FACE_CENTROIDS[fi])
        if d > best
            best = d
            best_fi = fi
        end
    end
    return best_fi
end

"""Return the 3 vertex indices (1-indexed) of face `fi`."""
_face_vertex_indices(fi::Int) = DGG_FACE_VERTEX_INDICES[fi]

#-----------------------------------------------------------------------------# Face Lambert Projection
# Each face has a Lambert azimuthal equal-area projection centered on its centroid.
# Local tangent-plane basis (u, v) oriented toward the first vertex of the face.

function _compute_face_bases()
    ico = icosahedron(:isea)
    verts_3d = [SVector{3,Float64}(v...) for v in ico.position]

    bases_u = Vector{SVector{3,Float64}}(undef, 20)
    bases_v = Vector{SVector{3,Float64}}(undef, 20)
    for fi in 1:20
        n = DGG_FACE_CENTROIDS[fi]
        # Orient u toward the first vertex of the face
        v1 = verts_3d[DGG_FACE_VERTEX_INDICES[fi][1]]
        t = v1 - dot(v1, n) * n  # project onto tangent plane
        u = normalize(t)
        v = cross(n, u)
        bases_u[fi] = u
        bases_v[fi] = v
    end
    return SVector{20}(bases_u), SVector{20}(bases_v)
end

const DGG_FACE_BASIS_U, DGG_FACE_BASIS_V = _compute_face_bases()

"""Project 3D unit vector `p` onto face `fi`'s Lambert tangent plane.
Returns (x, y) in the face's local coordinate system."""
function _face_lambert_forward(p::SVector{3,Float64}, fi::Int)
    n = DGG_FACE_CENTROIDS[fi]
    cos_a = clamp(dot(p, n), -1.0, 1.0)
    cos_a <= -1.0 + 1e-14 && return (0.0, 0.0)  # antipodal
    k = sqrt(2.0 / (1.0 + cos_a))
    x = k * dot(p, DGG_FACE_BASIS_U[fi])
    y = k * dot(p, DGG_FACE_BASIS_V[fi])
    return (x, y)
end

"""Inverse face Lambert: (x, y) in face `fi`'s tangent plane -> 3D unit vector."""
function _face_lambert_inverse(x::Float64, y::Float64, fi::Int)
    rho2 = x * x + y * y
    rho2 < 1e-30 && return DGG_FACE_CENTROIDS[fi]
    cos_a = 1.0 - rho2 / 2.0
    scale = sqrt(max(0.0, 1.0 - rho2 / 4.0))
    p = cos_a * DGG_FACE_CENTROIDS[fi] + scale * (x * DGG_FACE_BASIS_U[fi] + y * DGG_FACE_BASIS_V[fi])
    return normalize(p)
end
