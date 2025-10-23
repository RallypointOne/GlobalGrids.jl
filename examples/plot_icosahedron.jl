using GlobalGrids, GLMakie, GeoMakie

import GlobalGrids as GG
import GeometryBasics as GB


# Shapes
ico = GG.UnitIcosahedron{Float32}()

sph = GB.Sphere(ico)  # circumscribed sphere

# Plot
fig = Figure()
ax = Axis3(fig[1, 1], aspect = (1, 1, 1), title="Icosahedron and Circumscribed Sphere")
mesh!(ax, ico.mesh, color = :lightblue, shading = true)
wireframe!(ax, ico.mesh, color = :black, linewidth = 1)
mesh!(ax, sph, color = (:orange, 0.3), shading = true)

fig
