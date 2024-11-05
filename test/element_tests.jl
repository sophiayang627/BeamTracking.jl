using BeamTracking
using AcceleratorLattice
using Test

# define beams
e_minus = Species("electron")
mec2 = mass(e_minus)
ekin1 = 5.e3;  # eV
ekin2 = 1.e6;  # eV
ekin3 = 1.e9;  # eV
bg1 = sqrt(ekin1 / mec2 * (ekin1 / mec2 + 2))
bg2 = sqrt(ekin2 / mec2 * (ekin2 / mec2 + 2))
bg3 = sqrt(ekin3 / mec2 * (ekin3 / mec2 + 2))
pc1 = mec2 * bg1
pc2 = mec2 * bg2
pc3 = mec2 * bg3
zi1 = Beam(8; beta_gamma_0 = bg1)

# define elements
@ele start1 = BeginningEle(pc_ref = pc1, species_ref = e_minus);
@ele d1 = Drift(L = 0.15);
@ele start2 = BeginningEle(pc_ref = pc2, species_ref = e_minus);
@ele d2 = Drift(L = 0.75);
@ele start3 = BeginningEle(pc_ref = pc3, species_ref = e_minus);
@ele d3 = Drift(L = 2.00);

# define beamlines
bl1_drift = BeamLine([ start1, d1 ]);
bl2_drift = BeamLine([ start2, d2 ]);
bl3_drift = BeamLine([ start3, d3 ]);

# expand beamlines
lat1_drift = Lat([ bl1_drift ])
lat2_drift = Lat([ bl2_drift ])
lat3_drift = Lat([ bl3_drift ])

# test individual elements
@testset "element_tests" begin
  # drift
  @test lat1_drift["d1"][1].L == 0.15
  #@test track!(lat_drift["d1"][1], zf, zi);
  #      zf == [ x_f, px_f, y_f, py_f, z_f, pz_f ]

  # quadrupole

end

