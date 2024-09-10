using BeamTracking
using Test

@testset "BeamTracking.jl" begin
    beam0 = Beam(10000, d_x=Normal(0, 1e-3), d_px=Normal(0, 1e-4), d_pz=Normal(0, 1e-2))
    beamf = Beam(10000)
end
