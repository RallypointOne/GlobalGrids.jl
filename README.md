# GlobalGrids

[![Build Status](https://github.com/joshday/GlobalGrids.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/joshday/GlobalGrids.jl/actions/workflows/CI.yml?query=branch%3Amain)


## Technical Details on the Grid

- Represent earth with an icosahedron (12 vertices with 20 triangular faces).
- There are 12 base cells (pentagons) centered on the 12 vertices.

### Finding ZIndex from Lat/Lon

1. Convert `coord` to ECEF
2. Check `dot(face_enter, coord)`.  Maximum is the face you're in.
3. Rotate by face-specific rotaiton matrix.
4. 
