
using BeamTracking
using AcceleratorLattice
using Test

# # define elements
# @ele start = BeginningEle(pc_ref = 10.e6, species_ref = Species("electron"));
# @ele d1 = Drift(L = 0.55);
# @ele q1 = Quadrupole(L = 0.20 , Kn1 = -1.2)

# # define beamlines
# bl_drift = BeamLine("bl_drift", [ start, d1 ]);
# bl_quad = Beamline("bl_quad", [ start, q1 ])

# # expand beamlines
# lat_drift = expand([ bl_drift ])
# lat_quad = expand([ bl_quad ])



@testset "Linear" begin
#drift 
drift_ele = Linear.Drift(L=0.55)
#Quadrupole
quad_ele = Linear.Quadrupole(L=0.20, B1 = 1.0)

#single particle beam 
beami = Beam(px=1.0,py=1.0)
beamf = Beam()


Linear.track!(drift_ele, beamf, beami)
# Linear.track!(quad_ele, beamf, beami)

@test !(beamf == beami)
println(beamf,"linear")
# println(beamf, "quad")

end