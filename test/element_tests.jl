using BeamTracking
#using AtomicAndPhysicalConstants
using Test

# define elements
# -- drifts
l1 = 0.15
l2 = 0.75
l3 = 2.00
dr1 = MatrixKick.Drift(L=l1); #kwarg ctors
dr2 = MatrixKick.Drift(L=l2);
dr3 = MatrixKick.Drift(L=l3);
# -- quadrupoles


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
xi1  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200,  0.00200,  0.00200];
pxi1 = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075,  0.00075,  0.00075];
yi1  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100,  0.00100,  0.00100];
pyi1 = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030,  0.00030,  0.00030];
zi1  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000,  0.00000,  0.00000];
pzi1 = [ 0.001,  0.000, -0.001,  0.001,  0.00100,  0.00100,  0.00100,  0.00100];


zf1 = Coords(8);
zf1.x  .= [0.,                   0.,  0.,                  2.e-3,                2.112387648980866-3,  -2.112387648980866-3,  2.2247752979617317-3, 2.2247752979617317-3]
zf1.px .= [0.,                   0.,  0.,                  0.,                   7.5e-4,               -7.5e-4,               7.5e-4,               7.5e-4              ]
zf1.y  .= [0.,                   0.,  0.,                  1.e-3,                1.0449550595923462-3, -1.0449550595923462-3, 1.0899101191846924-3, 1.0899101191846926-3]
zf1.py .= [0.,                   0.,  0.,                  0.,                   3.e-4,                -3.e-4,                3.e-4,                3.e-4               ]
zf1.z  .= [1.4999999141696263-4, 0., -1.499999914426978-4, 1.4999999141696263-4, 1.4995115162148818-4,  1.4995115162148818-4, 2.9990230324297636-4, 2.9990230324297636-4]
zf1.pz .= [1.e-3,                0., -1.e-3,               1.e-3,                1.e-3,                 1.e-3,                1.e-3,                1.e-3               ]
beami1 = Bunch(e_minus, bg1, zi1)
beamf1 = Bunch(e_minus, bg1, zf1)


# test individual elements
@testset "element_tests" begin
  # drift
  @test dr1.L == l1
  @test dr2.L == l2
  #@test track!(beamf1, dr1, beami1);
  #      beamf.z == [ x_f, px_f, y_f, py_f, z_f, pz_f ]

  # quadrupole
#  @test norm(zf1.z - beamf1.v.z) < tol


end

