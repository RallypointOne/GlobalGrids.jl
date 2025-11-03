using GlobalGrids, GLMakie, GeoMakie, LinearAlgebra

import GlobalGrids as GG
import GeometryBasics as GB


#-----------------------------------------------------------------------------# plot_coordinates
function plot_coordinates()
    x = -180:5.0:180
	y = -90:5.0:90
	points = vec([GG.LonLat(Point2d((x,y))) for y in y, x in x])
    fig = Figure(size=(1000, 1000))

    grid = GG.Grid()
    verts = GG.LonLat.(GG.ECEF.(GG.R .* grid.ico.position))
    normals = GG.LonLat.(GG.ECEF.(GG.R .* grid.dual.position))

    poles = [GG.LonLat(Point2d(0.0, 90.0)), GG.LonLat(Point2d(0.0, -90.0))]
    equator = [GG.LonLat(Point2d(lon, 0.0)) for lon in -180:20:180]

    i = 0
    function f!(T::Type{<:GG.AbstractCoordinates{N}}) where {N}
        ax = N == 2 ?
            Axis(fig[i += 1, 1], title=string(T)) :
            Axis3(fig[i += 1, 1], title=string(T), aspect = (1, 1, 1))
        scatter!(ax, T.(points), color = :blue, markersize = 2)
        scatter!(ax, T.(verts), color = (:red, .5), markersize = 10, marker = :star8)
        # scatter!(ax, T.(normals), color = (:green, .5), markersize = 10, marker=:rect)
        text!(ax, T.(poles), text=["North", "South"], fontsize=10)
        scatter!(ax, T.(equator), color=:black, markersize=5)
        text!(ax, T.(equator), fontsize=10, text=string.(map(x -> x.pt[1], equator)))
    end

    f!(GG.LonLat)
    f!(GG.ECEF)
    f!(GG.ISEA)
    f!(GG.ISEACube)




    fig
end

#-----------------------------------------------------------------------------# plot_projections
function plot_base_cells(grid = Grid())
    fig = Figure()
    ax = Axis3(fig[1, 1], aspect = (1, 1, 1), title="Base Cell (Resolution = 0)")
    wireframe!(ax, GG.dual(grid.ico), color = :red, linewidth = 1)
    fig
end

#-----------------------------------------------------------------------------# icosahedron plots
# With circumscribed sphere
function plot_icosahedron()
    ico = GG.unit_icosahedron()
    sph = GG.unit_sphere()

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
        c = GG.centroid(tri)
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
