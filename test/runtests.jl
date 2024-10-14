using BeamTracking
using Distributions
using Test

@testset "BeamTracking.jl" begin
    beam0 = Beam(10000, d_x=Normal(0, 1e-3), d_px=Normal(0, 1e-4), d_pz=Normal(0, 1e-2))
    beamf = Beam(10000)

  #=
  @testset "element_tests" begin
    include("element_tests.jl")
  end
  =#
end

include("linear.jl")
