module GlobalGridsMakieExt

using GlobalGrids, Makie, GeoInterface

GeoInterface.@enable_makie Makie GlobalGrids.AbstractCoordinates
GeoInterface.@enable_makie Makie GlobalGrids.H3Cell

end
