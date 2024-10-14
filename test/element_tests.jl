using BeamTracking
using AcceleratorLattice
using Test

# define elements
@ele start = BeginningEle(pc_ref = 10.e6, species_ref = Species("electron"));
@ele d1 = Drift(L = 2.);

# define beamlines
bl_drift = BeamLine("bl_drift", [ start, d1 ]);

# expand beamlines
lat_drift = expand([ bl_drift ])

# test individual elements
@testset "element_tests" begin
  # drift
  @test lat_drift["d1"][1].L == 2.0
  #@test track!(lat_drift["d1"][1].L, zf, zi);
  #      zf == [ x_f, px_f, y_f, py_f, z_f, pz_f ]

  # quadrupole

end

