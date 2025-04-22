using BeamTracking, Beamlines, GTPSA, BenchmarkTools

function make_fodo(K1=0.40, L_quad=0.5, L_drift=5.0)
  qf = Quadrupole(K1=K1, L=L_quad, tracking_method=Linear())
  d1 = Drift(L=L_drift, tracking_method=Linear())
  qd = Quadrupole(K1=-qf.K1, L=L_quad, tracking_method=Linear())
  d2 = Drift(L=L_drift, tracking_method=Linear())
  return [qf, d1, qd, d2]
end

K1 = 0.40
L_quad = 0.5
L_drift = 5.0
N_fodo = 100

bl = Beamline([ele for i in 1:N_fodo for ele in make_fodo(K1,L_quad,L_drift)]; Brho_ref=60.0)
N_particle = 100

b0 = Bunch(N_particle, mem=BeamTracking.SoA)

# Track the bunch
track!(b0, bl)

# Also can do track! on individual elements
track!(b0, bl.line[1])

# And if you want to move the particle for-loop to the outside:
track!(b0, bl; outer_particle_loop=true)

# Can also track the bits representation:
bbl = BitsBeamline(bl)
track!(b0, bbl)

# GTPSA map:
const D = Descriptor(6,1)
v = @vars(D)
b0 = Bunch(v, mem=BeamTracking.AoS)

track!(b0, bl)
track!(b0, bbl)