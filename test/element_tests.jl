using BeamTracking
using AcceleratorLattice
using Test

# define elements
@ele start = BeginningEle(pc_ref = 10.e6, species_ref = Species("electron"));
@ele d1 = Drift(L = 2.);

# define beamlines
bl_drift = beamline("bl_drift", [ start, d1 ]);

# expand beamlines
lat_drift = expand("lat_drift", [ bl_drift ])

# test individual elements
@testset "element_tests" begin
  # drift
  @test lat_drift["d1"][1].L == 2.0
  #@test track!(d1, zf, zi) == [xz, pxf, yf, pyf, zf, pzf]

  # quadrupole

end

