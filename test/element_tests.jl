using BeamTracking
#using AtomicAndPhysicalConstants
using Test

# define elements
# -- drifts
l1 = 0.15
l2 = 0.75
l3 = 2.00
dr1 = MatrixKick.Drift(L = l1)
dr2 = MatrixKick.Drift(L = l2)
dr3 = MatrixKick.Drift(L = l3)
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
x1  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ]  #,  0.00200,  0.00200];
px1 = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ]  #,  0.00075,  0.00075];
y1  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ]  #,  0.00100,  0.00100];
py1 = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ]  #,  0.00030,  0.00030];
z1  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ]  #,  0.00000,  0.00000];
pz1 = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ]  #,  0.00100,  0.00100];

xf1  = [ 0.,  0.,                     0.,                   2.e-3,                 2.112387648980866e-3,  -2.112387648980866e-3  ]  #,  2.2247752979617317-3, 2.2247752979617317-3]
pxf1 = [ 0.,  0.,                     0.,                   0.,                    7.5e-4,                -7.5e-4                ]  #   7.5e-4,               7.5e-4              ]
yf1  = [ 0.,  0.,                     0.,                   1.e-3,                 1.0449550595923462e-3, -1.0449550595923462e-3 ]  #, 1.0899101191846924-3, 1.0899101191846926-3]
pyf1 = [ 0.,  0.,                     0.,                   0.,                    3.e-4,                 -3.e-4                 ]  #  3.e-4,                3.e-4               ]
zf1  = [ 0.,  1.4999999141696263e-4, -1.499999914426978e-4, 1.4999999141696263e-4, 1.4995115162148818e-4,  1.4995115162148818e-4 ]  #, 2.9990230324297636-4, 2.9990230324297636-4]
pzf1 = [ 0.,  1.e-3,                 -1.e-3,                1.e-3,                 1.e-3,                  1.e-3                 ]  #  1.e-3,                1.e-3               ]
beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
             x = x1, px = px1, y = y1, py = py1, z = z1, pz = pz1)


# test individual elements
@testset "element_tests" begin
  # drift
  @test dr1.L == l1
  @test dr2.L == l2
  track!(beam1, dr1);
  @test beam1.v.x  == xf1 #and beam1.px == pxf1 and beam1.y  == yf1 and beam1.py == pyf1 and beam1.z  == zf1 and beam1.pz == pzf1 

  # quadrupole
#  @test norm(zf1.z - beamf1.v.z) < tol


end

