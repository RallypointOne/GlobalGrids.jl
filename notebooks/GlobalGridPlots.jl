### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# ╔═╡ 756be6be-b1b7-11f0-b496-e997cfad5288
begin 
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	using GLMakie, GeoMakie, LinearAlgebra, CoordinateTransformations, Rotations, StaticArraysCore
	import GeoInterface as GI 
	import GeoFormatTypes as GFT 
	using GeometryBasics 
	import GeometryBasics as GB
	using GeometryBasics: Mesh
	using PlutoUI
	PlutoUI.TableOfContents()
end

# ╔═╡ 76bb5d17-2811-48fa-b9fa-37bbfb252e4f
md"# Functions"

# ╔═╡ 1ab5be94-9e8f-4645-b104-5df433fe744b
unit_sphere(T = Float64) = GB.Sphere(Point3{T}(0.0, 0.0, 0.0), 0.0)

# ╔═╡ 924f65b9-1320-41a3-90ff-3c176790d047
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

# ╔═╡ cf3d9aff-d702-40bd-8fed-2ff97327722d
centroid(x::GB.Ngon{D,T,N}) where {D,T,N} = sum(x.points) / N

# ╔═╡ 27c929d1-5bd0-411b-b397-a82ed628bb61
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

# ╔═╡ 20a78e55-33a1-45ae-bc6d-26d72b9fccea
function rotation_matrix(axis::SVector{3,T}, angle) where {T}
    K = SMatrix{3,3,T}(0, axis[3], -axis[2], -axis[3], 0, axis[1], axis[2], -axis[1], 0)
    return I + sin(angle) * K + (1 - cos(angle)) * (K * K)
end

# ╔═╡ 989834e6-c1c7-47cb-91d6-73b338f06113
md"# Icosahedron"

# ╔═╡ 6d27e665-9686-4fa7-8693-64dc7c02956f
ico = let
  φ = (1 + sqrt(5)) / 2
    verts = map(normalize, Point3d[
        (-1,  φ, 0),  (1,  φ, 0),  (-1, -φ, 0),  (1, -φ, 0),
        (0, -1,  φ), (0,  1,  φ), (0, -1, -φ), (0,  1, -φ),
        (φ,  0, -1), (φ,  0,  1), (-φ, 0, -1), (-φ, 0,  1)
    ])

    faces = GB.TriangleFace.([(1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12),  
		(2,6,10), (6,12,5), (12,11,3), (11,8,7), (8,2,9), (4,10,5), (4,5,3), 
		(4,3,7), (4,7,9), (4,9,10), (5,10,6), (3,5,12), (7,3,11), (9,7,8), (10,9,2)
	])
	GB.Mesh(verts, faces)
end

# ╔═╡ 8504b589-82ce-468c-af56-c1c3efe70d8e
function plot_polyhedron(ico = ico)
	fig = Figure()
	ax = Axis3(fig[1,1], aspect=(1,1,1))
	
	colors = distinguishable_colors(20, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    for (i, face) in enumerate(ico.faces)
        tri = ico[i]
        mesh!(ax, tri, color = colors[i], shading = NoShading)
		c = centroid(tri)
        text!(ax, 1.05GB.Point3f(c...), text = string(i), fontsize = 12,
              align = (:center, :center), color = :black)
    end

		
	wireframe!(ax, ico, color=:black)
	# wireframe!(ax, rotate(ico, RotXYZ(pi/6, pi / 6, pi / 6)))
	fig
end

# ╔═╡ ff15c322-e969-462e-8ffa-d87cd0c5f518
function plot_polyhedron!(ax, ico = ico)
	colors = distinguishable_colors(20, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    for (i, face) in enumerate(ico.faces)
        tri = ico[i]
        mesh!(ax, tri, color = colors[i], shading = NoShading)
		c = centroid(tri)
        text!(ax, 1.05GB.Point3f(c...), text = string(i), fontsize = 12,
              align = (:center, :center), color = :black)
    end
	
	wireframe!(ax, ico, color=:black)
	# wireframe!(ax, rotate(ico, RotXYZ(pi/6, pi / 6, pi / 6)))
end

# ╔═╡ a81f1beb-333d-4a3c-a70b-192ff9c03c5a
plot_polyhedron(ico)

# ╔═╡ 292e0c5a-d0f5-4630-bb72-26017b7cb10b
let 
    sph = unit_sphere()

    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = (1, 1, 1), title="Icosahedron and Circumscribed Sphere")
    mesh!(ax, ico, color = :lightblue, shading = true)
    wireframe!(ax, ico, color = :black, linewidth = 1)
    mesh!(ax, sph, color = (:orange, 0.3), shading = true)

    ax2 = Axis3(fig[1,2], aspect = (1, 1, 1), title="Icosahedron Projected onto Sphere")

    # The icosahedron vertices are already normalized (on the unit sphere)
    # Create a sphere mesh and overlay the icosahedron wireframe
    colors = distinguishable_colors(20, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    # Plot 3D icosahedron with colored faces
    for (i, face) in enumerate(ico.faces)
        tri = ico[i]
        mesh!(ax, tri, color = colors[i], shading = NoShading)
    end

    mesh!(ax2, sph, color = (:lightgray, 0.4), shading = true)

    # Add face numbers
    for (i, tri) in enumerate(ico)
        c = centroid(tri)
        text!(ax, 1.05GB.Point3f(c...), text = string(i), fontsize = 12,
              align = (:center, :center), color = :black)
    end

    # Function to interpolate along great circle arc between two points on unit sphere
    function great_circle_arc(p1, p2, n=50)
        # Angle between the two points
        θ = acos(clamp(dot(p1, p2), -1, 1))
        # Interpolate along the arc using slerp (spherical linear interpolation)
        [normalize((sin((1-t)*θ)/sin(θ)) * p1 + (sin(t*θ)/sin(θ)) * p2) for t in range(0, 1, length=n)]
    end

    # Plot the icosahedron edges as geodesic arcs on the sphere surface
    vertices = ico.position
    for face in ico.faces
        # Each face has 3 edges
        for i in 1:3
            v1 = vertices[face[i]]
            v2 = vertices[face[mod1(i+1, 3)]]
            arc = great_circle_arc(v1, v2, 50)
            lines!(ax2, arc, color = :blue, linewidth = 2)
        end
    end

    # Plot the vertices as points
    scatter!(ax2, vertices, color = :black, markersize = 10)
    # hidedecorations!(ax)
    # hidedecorations!(ax2)

    fig
end

# ╔═╡ a1f1785e-1fc0-40d1-96e4-35502f438639
md"## Dual"

# ╔═╡ 61ed0ee1-58b1-44dc-9856-9aae860d0b58
dual = let 
	m = ico
	verts = Point.(values(GB.face_normals(m.position, m.faces)))
	face_normals = m.position
	faces = GB.NgonFace{5,Int}[]
	map(face_normals) do n
		# Find the 5 vertices closest to this normal
		pentagon = sort(verts, by = v -> -dot(v, n))[1:5]
		# Make sure vertices are counter-clockwise
		counterclockwise!(pentagon)

		idx = [findfirst(==(v), verts) for v in pentagon]
		push!(faces, GB.NgonFace(idx...))
	end
	Mesh(verts, faces)
end

# ╔═╡ 0fe3aaa4-1ab7-487f-af06-ac71e2513459
plot_polyhedron(dual)

# ╔═╡ b4ee7aad-4d89-4ab7-bd35-7da5abcb58f2
begin
	function rotate(m::Mesh, R::Rotations.Rotation{3})
	    new_verts = normalize.([Point(Tuple(R * GB.Vec(p))) for p in coordinates(m)])
	    GB.Mesh(new_verts, m.faces)
	end
	
	rotate(m::Mesh, ::typeof(I)) = m
	
	function rotate(m::Mesh; x=0, y=0, z=0)
		Rx = iszero(x) ? I : RotX(x)
		Ry = iszero(y) ? I : RotY(x) 
		Rz = iszero(z) ? I : RotZ(x)
		R = Rx * Ry * Rz 
		rotate(m, R)
	end
end

# ╔═╡ c3ff21dc-fb9a-4fb1-b598-042b9cf94659
md"# Rotations"

# ╔═╡ 8b0324be-a958-45cd-885d-cd928ec375b5
function animate_rotate(mesh, axis=:x; framerate=30)
	θ = Observable(0.0)

	fig = Figure()
	θstring = @lift string(round(rad2deg($θ), digits=0), '°')
	ax = Axis3(fig[1,1], aspect=(1,1,1), title=@lift("Rotation around $axis: $($θstring)"))
	Rot = if axis == :x 
		@lift RotX($θ)	
	elseif axis == :y 
		@lift RotY($θ)
	elseif axis == :z 
		@lift RotZ($θ)
	else
		error("Axis must be :x, :y, or :z")
	end
	points = @lift rotate(mesh, $Rot).position
	poly = @lift Mesh($points, mesh.faces)

	wireframe!(ax, poly, color=:black)
	mesh!(ax, poly, shading = NoShading, color=1:12)

	normals = @lift 1.01Point.(values(GB.face_normals($poly.position, $poly.faces)))
	text!(ax, normals, align=(:center, :center), text=string.(1:20))

	θs = 1:2:359
	file = record(fig, "animate_$axis.mp4", θs; framerate) do angle 
		θ[] = deg2rad(angle)
	end
	PlutoUI.LocalResource(file)
end

# ╔═╡ 4bded207-9432-435d-b0ce-8e4c640b0db4
animate_rotate(rotate(ico, RotY(pi/3)), :z)

# ╔═╡ 34f7ae24-89a7-472e-bb30-a4d81c419c82
md"# ISEA"

# ╔═╡ 1ab0f673-8575-4539-9ba4-1d3f60229917
md"## Orientations"

# ╔═╡ e13c7bae-4d9f-40f4-904d-f7d6b1fc7013
ORIENTATION_POLES = (vertex = Point3d(0,0, 1), azimuth = 0.0)

# ╔═╡ e005131e-580c-468b-8206-928178cd7a56
# vertex = normalize(GG.ECEF(GG.LonLat(11.25, 58.282525590)).pt)
ORIENTATION_ISEA = (vertex = Point3d(0.518136909297074, 0.1030638392550601, 0.8490653615959627), azimuth = 0.0)

# ╔═╡ 95508cff-c8cb-46b2-b74c-e70ea555c5ed
# vertex = normalize(GG.ECEF(GG.LonLat(-5.245390, 2.3008820)).pt)
ORIENTATION_DYMAXION = (vertex = Point3d(0.9950201403480873, -0.09134877248651022, 0.039878842346294276), azimuth = deg2rad(7.466580))

# ╔═╡ 5ff781eb-218a-4eec-b0b7-0356d06df496
md"## Create Icosahedrons per Orientations"

# ╔═╡ 6aeaa0e5-3331-4a53-9134-db89c890e739
function isea_icosahedron(; vertex::Point3{T} = Point3d(0,0,1), azimuth::T = 0.0) where {T}
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

	R = 6_371_007.180918475

    new_vertices = R .* normalize.([Point{3,T}((Rot * SVector{3,T}(v...))...) for v in unit_ico.position])

    return Mesh(new_vertices, unit_ico.faces)
end

# ╔═╡ 2221d96f-c0bb-476c-85e9-b3f7ca7572df
ico_poles = isea_icosahedron(; ORIENTATION_POLES...)

# ╔═╡ dbe42cc6-cd12-4478-a176-a23881fca3ee
ico_isea = isea_icosahedron(; ORIENTATION_ISEA...)

# ╔═╡ 3c1d0fd4-a7b3-4068-aea4-0a8b33d2e64a
ico_dymaxion = isea_icosahedron(; ORIENTATION_DYMAXION...)

# ╔═╡ 26b7089f-049f-40f4-801b-d8bce6785e49


# ╔═╡ Cell order:
# ╟─756be6be-b1b7-11f0-b496-e997cfad5288
# ╟─76bb5d17-2811-48fa-b9fa-37bbfb252e4f
# ╟─1ab5be94-9e8f-4645-b104-5df433fe744b
# ╟─924f65b9-1320-41a3-90ff-3c176790d047
# ╠═cf3d9aff-d702-40bd-8fed-2ff97327722d
# ╟─27c929d1-5bd0-411b-b397-a82ed628bb61
# ╟─8504b589-82ce-468c-af56-c1c3efe70d8e
# ╟─ff15c322-e969-462e-8ffa-d87cd0c5f518
# ╟─20a78e55-33a1-45ae-bc6d-26d72b9fccea
# ╟─989834e6-c1c7-47cb-91d6-73b338f06113
# ╟─6d27e665-9686-4fa7-8693-64dc7c02956f
# ╟─a81f1beb-333d-4a3c-a70b-192ff9c03c5a
# ╟─292e0c5a-d0f5-4630-bb72-26017b7cb10b
# ╟─a1f1785e-1fc0-40d1-96e4-35502f438639
# ╟─61ed0ee1-58b1-44dc-9856-9aae860d0b58
# ╟─0fe3aaa4-1ab7-487f-af06-ac71e2513459
# ╟─b4ee7aad-4d89-4ab7-bd35-7da5abcb58f2
# ╟─c3ff21dc-fb9a-4fb1-b598-042b9cf94659
# ╟─8b0324be-a958-45cd-885d-cd928ec375b5
# ╠═4bded207-9432-435d-b0ce-8e4c640b0db4
# ╟─34f7ae24-89a7-472e-bb30-a4d81c419c82
# ╟─1ab0f673-8575-4539-9ba4-1d3f60229917
# ╟─e13c7bae-4d9f-40f4-904d-f7d6b1fc7013
# ╟─e005131e-580c-468b-8206-928178cd7a56
# ╟─95508cff-c8cb-46b2-b74c-e70ea555c5ed
# ╟─5ff781eb-218a-4eec-b0b7-0356d06df496
# ╟─6aeaa0e5-3331-4a53-9134-db89c890e739
# ╟─2221d96f-c0bb-476c-85e9-b3f7ca7572df
# ╟─dbe42cc6-cd12-4478-a176-a23881fca3ee
# ╟─3c1d0fd4-a7b3-4068-aea4-0a8b33d2e64a
# ╠═26b7089f-049f-40f4-801b-d8bce6785e49
