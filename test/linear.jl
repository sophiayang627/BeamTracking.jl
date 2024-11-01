
using BeamTracking
using Test
using GTPSA


@testset "Linear" begin
#drift 
drift_ele = Linear.Drift(L=0.55)
#Quadrupole
quad_ele1 = Linear.Quadrupole(L=0.55, B1 = 0.1) #positive B1
quad_ele2 = Linear.Quadrupole(L=0.55, B1 = 0.0) #no strength, same as drift
quad_ele3 = Linear.Quadrupole(L=0.55, B1 = -0.1) #negative B1
quad_ele4 = Linear.Quadrupole(L=0.55, B1 = 2.0) #big B1

#single particle beam
d = Descriptor(6,1)
bi = Beam(x=1.0,px=0.1,y=1.0,py=0.1,z=0.1,pz=0.1)
bf = Beam()


#Drift 
Linear.track!(drift_ele, bf, bi)
@test !(bf == bi)


println(bf,"drift")
par = Particle(b=bf)
mat6 = GTPSA.jacobian(par.z)
println(bf,mat6,"drift")

 # if mat6 = mat from bmad

#Quad
Linear.track!(quad_ele2, bf, bi) # positive B1
@test !(bf == bi)
println(bf,"quad")
#par = Particle(b=bf)
#mat6=GTPSA.jacobian(par.z)
#println(bf,mat6,"quad with k1=0")

#check horizontal focusing, vertical defocusing for B>0
#check vertical focusing, horizontal defocusing for B<0

#0 strength quadrupole should give the same result as a drift
end