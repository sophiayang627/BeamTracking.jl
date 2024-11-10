using Test,
      BeamTracking,
      Distributions,
      JET,
      GTPSA

@testset "BeamTracking.jl" begin
    # Linear drift test:
    tol = 1e-14
    d = Descriptor(6, 1) # 6 variables, order 1
    beam = Beam(d, beta_gamma_ref=1.0) # Creates Beam as identity map using GTPSA Descriptor
    gamma_ref = sr_gamma(beam.beta_gamma_ref)
    L_d = 5.0
    drift = Linear.Drift(L=L_d)
    track!(beam, drift)

    p = Particle(beam, 1) # Converts SoA to 1 struct
    
    M_drift_expected = [1.0  L_d  0.0  0.0  0.0  0.0;
                        0.0  1.0  0.0  0.0  0.0  0.0;
                        0.0  0.0  1.0  L_d  0.0  0.0;
                        0.0  0.0  0.0  1.0  0.0  0.0;
                        0.0  0.0  0.0  0.0  1.0  L_d/gamma_ref^2;
                        0.0  0.0  0.0  0.0  0.0  1.0]
    
    # 1) Correctness:
    @test norm(GTPSA.jacobian(p.v) - M_drift_expected) < tol
    # 2) Type-stability:
    @test_opt track!(beam, drift) 
    # 3) Allocations:
    @test @allocations(track!(beam, drift)) == 0
end

#include("linear.jl")
