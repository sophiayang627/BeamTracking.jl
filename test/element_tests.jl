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
xf_dr1  = [ 0.,  0.,                     0.,                     2.e-3,                  2.1123876489808660e-3, -2.1123876489808660e-3 ]
pxf_dr1 = [ 0.,  0.,                     0.,                     0.,                     7.5e-4,                -7.5e-4                ]
yf_dr1  = [ 0.,  0.,                     0.,                     1.e-3,                  1.0449550595923462e-3, -1.0449550595923462e-3 ]
pyf_dr1 = [ 0.,  0.,                     0.,                     0.,                     3.e-4,                 -3.e-4                 ]
zf_dr1  = [ 0.,  1.4710284465059844e-4, -1.4711135596840458e-4,  1.4710284465059844e-4,  1.4705400485512816e-4,  1.4705400485512816e-4 ]
pzf_dr1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qf1.final
xf_qf1  = [ 0.,  0.,                     0.,                     4.4834792779340600e-3,  4.5356143287504990e-3, -4.5356143287504990e-3 ]
pxf_qf1 = [ 0.,  0.,                     0.,                     1.1607821136240948e-1,  1.1776162121669208e-1, -1.1776162121669208e-1 ]
yf_qf1  = [ 0.,  0.,                     0.,                     1.2400905948673489e-4,  1.3427609296030678e-4, -1.3427609296030678e-4 ]
pyf_qf1 = [ 0.,  0.,                     0.,                    -2.8691666098954356e-2, -2.8653744321335432e-2,  2.8653744321335432e-2 ]
zf_qf1  = [ 0.,  4.903428155019947e-5,  -4.903711865613486e-5,  -4.8701323139842656e-5, -5.1605970700562340e-5, -5.1605970700562340e-5 ]
pzf_qf1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]

# beam1.qd1.final
xf_qd1  = [ 0.,  0.,                     0.,                     2.4834727340278294e-4,  2.7409827330240064e-4, -2.7409827330240064e-4 ]
pxf_qd1 = [ 0.,  0.,                     0.,                    -5.7391734255615580e-2, -5.7299059109537503e-2,  5.7299059109537503e-2 ]
yf_qd1  = [ 0.,  0.,                     0.,                     2.2414071638535435e-3,  2.2621899982980444e-3, -2.2621899982980444e-3 ]
pyf_qd1 = [ 0.,  0.,                     0.,                     5.8033153149417160e-2,  5.8705242601074584e-2, -5.8705242601074584e-2 ]
zf_qd1  = [ 0.,  4.9034281550199470e-5, -4.9037118656134860e-5, -1.1242623339116514e-5, -1.1118475705649755e-5, -1.1118475705649755e-5 ]
pzf_qd1 = [ 0.,  1.e-3,                 -1.e-3,                  1.e-3,                  1.e-3,                  1.e-3                 ]


# test individual elements
@testset "element_tests" begin
  # === drifts ===
  beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
               x = copy(x1), px = copy(px1), y = copy(y1), py = copy(py1), z = copy(z1), pz = copy(pz1))
  track!(beam1, dr1);
  @test beam1.v.x  == xf_dr1 && beam1.v.px == pxf_dr1 &&
        beam1.v.y  == yf_dr1 && beam1.v.py == pyf_dr1 &&
        beam1.v.z  == zf_dr1 && beam1.v.pz == pzf_dr1

  # === quadrupoles ===
  # 5 keV
  beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
               x = copy(x1), px = copy(px1), y = copy(y1), py = copy(py1), z = copy(z1), pz = copy(pz1))
  track!(beam1, qf1);
  @test beam1.v.x  ≈  xf_qf1  (rtol=5.e-13)
  @test beam1.v.px ≈  pxf_qf1 (rtol=5.e-13)
  @test beam1.v.y  ≈  yf_qf1  (rtol=5.e-13)
  @test beam1.v.py ≈  pyf_qf1 (rtol=5.e-13)
  @test beam1.v.z  ≈  zf_qf1  (rtol=5.e-13)
  @test beam1.v.pz == pzf_qf1
  beam1 = Beam(species = e_minus, beta_gamma_ref = bg1,
               x = copy(x1), px = copy(px1), y = copy(y1), py = copy(py1), z = copy(z1), pz = copy(pz1))
  track!(beam1, qd1);
  @test beam1.v.x  ≈  xf_qd1  (rtol=5.e-13)
  @test beam1.v.px ≈  pxf_qd1 (rtol=5.e-13)
  @test beam1.v.y  ≈  yf_qd1  (rtol=5.e-13)
  @test beam1.v.py ≈  pyf_qd1 (rtol=5.e-13)
  @test beam1.v.z  ≈  zf_qd1  (rtol=5.e-13)
  @test beam1.v.pz == pzf_qd1


end

