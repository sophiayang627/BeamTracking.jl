using BeamTracking
using AtomicAndPhysicalConstants
using Test

# define beams
# -- species
e_minus = Species("electron")
mec2 = 0.51099895069e6; #massof(e_minus)
# -- kinetic energy
ek1 = 5.e3;  # eV
ek2 = 1.e6;  # eV
ek3 = 1.e9;  # eV
# -- βγ
bg1 = sqrt(ek1 / mec2 * (ek1 / mec2 + 2))
bg2 = sqrt(ek2 / mec2 * (ek2 / mec2 + 2))
bg3 = sqrt(ek3 / mec2 * (ek3 / mec2 + 2))
# -- pc
pc1 = mec2 * bg1
pc2 = mec2 * bg2
pc3 = mec2 * bg3
# -- beams
beami1 = Beam(8; beta_gamma_0 = bg1)

# define elements
# -- drifts
dr1 = Symplectic.Drift(0.15);
dr2 = Symplectic.Drift(0.75);
dr3 = Symplectic.Drift(2.00);
# -- quadrupoles


# test individual elements
@testset "element_tests" begin
  # drift
  @test dr1.L == 0.15
  @test dr2.L == 0.75
  #@test track!(d1, beamf, beami);
  #      beamf.z == [ x_f, px_f, y_f, py_f, z_f, pz_f ]

  # quadrupole

end

