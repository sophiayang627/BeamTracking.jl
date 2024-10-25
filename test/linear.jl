
using BeamTracking
using Test


@testset "Linear" begin
#drift 
drift_ele = Linear.Drift(L=0.55)
#Quadrupole
quad_ele = Linear.Quadrupole(L=0.20, B1 = 1.0)

#single particle beam 
bi1 = Beam(x=1.0,px=1.0,py=1.0)
bf1 = Beam()

bi2 = Beam(x=1.0,px=1.0,py=1.0)
bf2 = Beam()


Linear.track!(drift_ele, bf1, bi1)

Linear.track!(quad_ele, bf2, bi2)

@test !(bf1 == bi1)
@test !(bf2 == bi2)

println(bf1,"linear")
println(bf2, "quad")

end