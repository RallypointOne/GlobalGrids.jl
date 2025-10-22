module GlobalGrids

import GeoInterface as GI
import GeometryOps as GO
import GeoFormatTypes as GFT
import GeometryBasics as GB
import StyledStrings: @styled_str

using LinearAlgebra


export IGEO7

#-----------------------------------------------------------------------------# utils
"Approximate radius of the earth in meters (WGS84)."
const R = 6_371_000



#-----------------------------------------------------------------------------# Icosahedron
struct Icosahedron{T}
    mesh::GB.Mesh{3, T, GB.TriangleFace{Int}}

    function Icosahedron{T}(radius=1.0) where {T}
        φ = (1 + sqrt(5)) / 2
        verts = GB.Point{3,T}[
            (-1,  φ, 0),  (1,  φ, 0),  (-1, -φ, 0),  (1, -φ, 0),
            (0, -1,  φ), (0,  1,  φ), (0, -1, -φ), (0,  1, -φ),
            (φ,  0, -1), (φ,  0,  1), (-φ, 0, -1), (-φ, 0,  1)
        ]
        verts = GB.Point{3,T}.(radius .* normalize.(verts))

        faces = GB.TriangleFace.([(1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12), (2,6,10), (6,12,5),
            (12,11,3), (11,8,7), (8,2,9),  (4,10,5),  (4,5,3), (4,3,7),  (4,7,9),  (4,9,10), (5,10,6),
            (3,5,12), (7,3,11), (9,7,8),  (10,9,2)
        ])
        return new{T}(GB.Mesh(verts, faces))
    end
end
Icosahedron(radius=1.0) = Icosahedron{Float64}(radius)

radius(ico::Icosahedron) = norm(ico.mesh.position[1])

function Base.show(io::IO, o::Icosahedron)
    print(io, styled"{bright_cyan:Icosahedron}(radius = $(radius(o)))")
end



include("IGEO7.jl")

end
