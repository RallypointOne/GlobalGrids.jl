module GeoConstants

# Copied from Snyder (1992) "An Equal-Area Map Projection for Polyhedral Globes":
#
# g = spherical distance, degrees, from center of polygon face to any of its vertices on the globe.
# G = spherical angle, degrees, between radius vector to center and adjacent edge of spherical polygon
# on the globe.
# θ = plane angle, degrees, between radius vector to center and adjacent edge of plane polygon.
# ω = maximum value of maximum angular deformation, degrees, occurring along a radius to each
# vertex, but at the center on all equal-area polygons except the hexagons of the truncated icosahedron,
# where it occurs at each vertex. On the Gnomonic projection, this occurs at each vertex.
# a = maximum value of the maximum scale factor, occurring where ω is maximum.
# b = minimum value of the minimum scale factor, occurring where ω is maximum on equal-area
# polygons and at the centers of all Gnomonic polygons.
# NOTE: The values of a and b on the pentagon faces of the equal-area truncated icosahedron reflect the
# higher area scale on those faces, taking the hexagon faces at true area scale. On the Gnomonic faces,
# all faces are considered tangent to the sphere except for the pentagons of the truncated icosahedron.
hexagon_constants = (;
    g = 37.37736814,
    G = 36.0,
    θ = 30.0,
    ω = 17.27,
    a = 1.163,
    b = 0.860
)


"Authalic radius of the earth in meters (WGS84)"
const R = 6_371_007.180918475

"Equitorial radius (semi-major axis) of the earth (WGS84)."
const A = 6_378_137

"Polar radius (semi-minor axis) of the earth (WGS84)."
const B = 6_356_752

"Flattening of the earth (WGS84)."
const F = 1 - (B // A)  # flattening

"First eccentricity squared of the earth (WGS84)."
const E2 = 1 - A ^ 2 // B ^ 2  # first eccentricity squared

"Second eccentricity squared of the earth (WGS84)."
const E′2 = A ^ 2 // B ^ 2 - 1  # second eccentricity squared

"""
    prime_vertical_radius(lat°)

Reference:  https://en.wikipedia.org/wiki/Earth_radius#Prime_vertical
"""
prime_vertical_radius(lat) = A ^ 2 / sqrt(A ^ 2 * cosd(lat) ^ 2 + B ^ 2 * sind(lat) ^ 2)

end
