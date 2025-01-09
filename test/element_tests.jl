using BeamTracking
#using AtomicAndPhysicalConstants
using Test

# define elements
# -- drifts
ld1 = 0.15;  # m
ld2 = 0.75;  # m
ld3 = 2.00;  # m
dr1 = MatrixKick.Drift(L = ld1)
dr2 = MatrixKick.Drift(L = ld2)
dr3 = MatrixKick.Drift(L = ld3)

# -- quadrupoles
lq1 = 0.05;  # m
lq2 = 0.30;  # m
lq3 = 0.75;  # m
gr1 = 0.20;  # T/m
qf1 = MatrixKick.Quadrupole(L = lq1, Bn1 =  gr1)
qd1 = MatrixKick.Quadrupole(L = lq1, Bn1 = -gr1)


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
# beam1.initial
x1  = [ 0.000,  0.000,  0.000,  0.002,  0.00200, -0.00200 ]  #,  0.00200,  0.00200];
px1 = [ 0.000,  0.000,  0.000,  0.000,  0.00075, -0.00075 ]  #,  0.00075,  0.00075];
y1  = [ 0.000,  0.000,  0.000,  0.001,  0.00100, -0.00100 ]  #,  0.00100,  0.00100];
py1 = [ 0.000,  0.000,  0.000,  0.000,  0.00030, -0.00030 ]  #,  0.00030,  0.00030];
z1  = [ 0.000,  0.000,  0.000,  0.000,  0.00000,  0.00000 ]  #,  0.00000,  0.00000];
pz1 = [ 0.000,  0.001, -0.001,  0.001,  0.00100,  0.00100 ]  #,  0.00100,  0.00100];

# beam1.dr1.final
xf1_dr1  = [ 0.,  0.,                     0.,                     2.e-3,                  2.1123876489808660e-3, -2.1123876489808660e-3 ]
pxf1_dr1 = [ 0.,  0.,                     0.,                     0.,                     7.5e-4,                -7.5e-4                ]
yf1_dr1  = [ 0.,  0.,                     0.,                     1.e-3,                  1.0449550595923462e-3, -1.0449550595923462e-3 ]
pyf1_dr1 = [ 0.,  0.,                     0.,                     0.,                     3.e-4,                 -3.e-4                 ]
zf1_dr1  = [ 0.,  1.4710284465059844e-4, -1.4711135596840458e-4,  1.4710284465059844e-4,  1.4705400485512816e-4,  1.4705400485512816e-4 ]
pzf1_dr1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qf1.final
xf1_qf1   = [0.,  0.,                     0.,                     2.4834727340278294e-4,  2.7409827330240064e-4, -2.7409827330240064e-4 ]
pxf1_qf1  = [0.,  0.,                     0.,                    -5.7391734255615580e-2, -5.7299059109537503e-2,  5.7299059109537503e-2 ]
yf1_qf1   = [0.,  0.,                     0.,                     2.2414071638535435e-3,  2.2621899982980444e-3, -2.2621899982980444e-3 ]
pyf1_qf1  = [0.,  0.,                     0.,                     5.8033153149417160e-2,  5.8705242601074584e-2, -5.8705242601074584e-2 ]
zf1_qf1   = [0.,  4.9034281550199470e-5, -4.9037118656134860e-5, -1.1242623339116514e-5,  1.1118475705649755e-5, -1.1118475705649755e-5 ]
pzf1_qf1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qd1.final

# test individual elements
@testset "element_tests" begin
  # drift
  beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
               x = x1, px = px1, y = y1, py = py1, z = z1, pz = pz1)
  track!(beam1, dr1);
  @test beam1.v.x  == xf1_dr1 && beam1.v.px == pxf1_dr1 &&
        beam1.v.y  == yf1_dr1 && beam1.v.py == pyf1_dr1 &&
        beam1.v.z  == zf1_dr1 && beam1.v.pz == pzf1_dr1

  # quadrupole
  beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
               x = x1, px = px1, y = y1, py = py1, z = z1, pz = pz1)
  track!(beam1, qf1);
  @test beam1.v.x  == xf1_qf1 && beam1.v.px == pxf1_qf1 &&
        beam1.v.y  == yf1_qf1 && beam1.v.py == pyf1_qf1 &&
        beam1.v.z  == zf1_qf1 && beam1.v.pz == pzf1_qf1
#  @test norm(zf1.z - beamf1.v.z) < tol


end

