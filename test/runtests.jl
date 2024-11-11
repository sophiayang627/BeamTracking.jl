using Test,
      BeamTracking,
      Distributions,
      JET,
      BenchmarkTools,
      GTPSA

@testset "BeamTracking.jl" begin
    # Test setup
    BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
    BenchmarkTools.DEFAULT_PARAMETERS.evals = 2

    # Linear tests -----------------------
    d = Descriptor(6, 1) # 6 variables, order 1
    beta_gamma_ref = 1.0
    gamma_ref = sr_gamma(beta_gamma_ref)
    species = Species("electron")
    brho_ref = brho(massof(species), beta_gamma_ref, chargeof(species)) 
    tol = 1e-14

    # Linear drift test:
    beam = Beam(d, beta_gamma_ref=beta_gamma_ref) # Creates Beam as identity map using GTPSA Descriptor
    L_d = 5.0
    drift = Linear.Drift(L=L_d)
    track!(beam, drift)
    
    M_drift_expected = [1.0  L_d  0.0  0.0  0.0  0.0;
                        0.0  1.0  0.0  0.0  0.0  0.0;
                        0.0  0.0  1.0  L_d  0.0  0.0;
                        0.0  0.0  0.0  1.0  0.0  0.0;
                        0.0  0.0  0.0  0.0  1.0  L_d/gamma_ref^2;
                        0.0  0.0  0.0  0.0  0.0  1.0]
    
    p = Particle(beam, 1) # Converts SoA to 1 struct

    # 1) Correctness:
    @test norm(GTPSA.jacobian(p.v) - M_drift_expected) < tol
    # 2) Type-stability:
    @test_opt track!(beam, drift) 
    # 3) Allocations:
    @test @ballocated(track!($beam, $drift)) == 0

    # Linear Quadrupole test:
    # Focusing:
    beam = Beam(d, beta_gamma_ref=beta_gamma_ref)
    L_q = 1.2
    K1n = 0.36
    B1 = K1n*brho_ref
    qf = Linear.Quadrupole(B1=B1,L=L_q)
    work = BeamTracking.get_work(beam, Val{1}())
    track!(beam, qf, work=work)

    p = Particle(beam, 1) # Converts SoA to 1 struct
    
    M_qf_x = [cos(sqrt(K1n)*L_q)            sincu(sqrt(K1n)*L_q)*L_q;  
              -sqrt(K1n)*sin(sqrt(K1n)*L_q) cos(sqrt(K1n)*L_q)      ;]
    M_qf_y = [cosh(sqrt(K1n)*L_q)           sinhcu(sqrt(K1n)*L_q)*L_q; 
              sqrt(K1n)*sinh(sqrt(K1n)*L_q) cosh(sqrt(K1n)*L_q)      ;]
    M_qf_z = [1.0 L_q/gamma_ref^2;
              0.0 1.0            ;]

    M_qf_expected = zeros(6,6)
    M_qf_expected[1:2,1:2] = M_qf_x
    M_qf_expected[3:4,3:4] = M_qf_y
    M_qf_expected[5:6,5:6] = M_qf_z
    # 1) Correctness:
    @test norm(GTPSA.jacobian(p.v) - M_qf_expected) < tol
    # 2) Type-stability:
    @test_opt track!(beam, qf, work=work) 
    # 3) Allocations:
    @test @ballocated(track!($beam, $qf, work=$work)) == 0

    # Defocusing:
    beam = Beam(d, beta_gamma_ref=beta_gamma_ref)
    qd = Linear.Quadrupole(B1=-B1,L=L_q)
    track!(beam, qd, work=work)

    p = Particle(beam, 1)

    M_qd_expected = zeros(6,6)
    M_qd_expected[1:2,1:2] = M_qf_y
    M_qd_expected[3:4,3:4] = M_qf_x
    M_qd_expected[5:6,5:6] = M_qf_z

    # 1) Correctness:
    @test norm(GTPSA.jacobian(p.v) - M_qd_expected) < tol
    # 2) Type-stability:
    @test_opt track!(beam, qd, work=work) 
    # 3) Allocations:
    @test @ballocated(track!($beam, $qd, work=$work)) == 0
end

#include("linear.jl")
