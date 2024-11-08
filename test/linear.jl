
using BeamTracking
using Test
using GTPSA


@testset "Linear" begin
#elements
 #drift
  drift_ele = Linear.Drift(L=0.50)
 #Quadrupole
  quad_ele1 = Linear.Quadrupole(L=0.4, B1 = 0.0020) #positive B1, quad defocusing for e-
  quad_ele2 = Linear.Quadrupole(L=0.4, B1 = -0.0020) #negative B1, quad focusing for e-
  quad_ele3 = Linear.Quadrupole(L=0.4, B1 = 2.0) #big B1
  quad_ele4 = Linear.Quadrupole(L=0.4, B1 = 0.0) #no strength, same as drift

 #Sbend
  #sb_ele1 = Linear.SBend(L = 1.5, B0 = 0.000446) # positive field

#single particle beam
 d = Descriptor(6,1)
 bi = Beam(d)
 bf_d = Beam(d)
 bf_q1 = Beam(d)
 bf_q2 = Beam(d)
 bf_q3 = Beam(d)
 bf_q4 = Beam(d)
 bf_s1 = Beam(d)

#Drift 
 Linear.track!(drift_ele, bf_d, bi)
 @test !(bf_d == bi)
 par_d = Particle(b=bf_d)
 mat_d = GTPSA.jacobian(par_d.z)
 println("\n","drift")
 show(stdout, "text/plain", mat_d)



#Quad Defocusing for electron
 Linear.track!(quad_ele1, bf_q1, bi) # Positive B1
 @test !(bf_q1 == bi)
 par_q1 = Particle(b=bf_q1)
 mat_q1 = GTPSA.jacobian(par_q1.z)
 println("\n","quad with k1=-1.2, Defocusing e-")
 show(stdout, "text/plain", mat_q1)

#Quad Focusing 
 Linear.track!(quad_ele2, bf_q2, bi) # negative B1
 @test !(bf_q2 == bi)
 par_q2 = Particle(b=bf_q2)
 mat_q2 = GTPSA.jacobian(par_q2.z)
 println("\n","quad with k1=1.2, Focusing e-")
 show(stdout, "text/plain", mat_q2)

#Large k, defocusing
 Linear.track!(quad_ele3, bf_q3, bi) # negative B1
 @test !(bf_q3 == bi)
 par_q3 = Particle(b=bf_q3)
 mat_q3 = GTPSA.jacobian(par_q3.z)
 println("\n","quad with B=2, defocusing e-")
 show(stdout, "text/plain", mat_q3)

#0 strength quadrupole should give the same result as a drift
 Linear.track!(quad_ele4, bf_q4, bi) # negative B1
 @test !(bf_q4 == bi)
 par_q4 = Particle(b=bf_q4)
 mat_q4 = GTPSA.jacobian(par_q4.z)
 println("\n","quad with B1=0, simplify to drift")
 show(stdout, "text/plain", mat_q4)

#Sbend
    # Linear.track!(sb_ele1, bf_s1, bi) # negative B1
    #  @test !(bf_s1 == bi)
    #  par_s1 = Particle(b=bf_s1)
    #  mat_s1 = GTPSA.jacobian(par_s1.z)
    #  println("\n","Sbend,positive B0")
    #  show(stdout, "text/plain", mat_s1)

end