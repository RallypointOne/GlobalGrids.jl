using GlobalGrids
using Test

import GlobalGrids as GG
import GeometryBasics as GB

@testset "GlobalGrids.jl" begin

@testset "Coordinates" begin
    lldeg = GG.LonLatDeg(GB.Point2(-75.0, 45.0))
    llrad = GG.LonLatRad(lldeg)
    isea = GG.ISEA(llrad)
    lldeg2 = GG.LonLatDeg(isea)

end # Coordinates
end
