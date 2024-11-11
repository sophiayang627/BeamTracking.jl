using Test,
      BeamTracking,
      Distributions,
      JET,
      BenchmarkTools,
      GTPSA

BenchmarkTools.DEFAULT_PARAMETERS.gctrial = false
BenchmarkTools.DEFAULT_PARAMETERS.evals = 2
  
# Soon we generalize this to test_map...
function test_matrix(ele, n_work, M_expected; type_stable=true, no_allocs=true, tol=1e-14, 
                                              beta_gamma_ref=1.0, species=Species("electron")   )
    
  d = Descriptor(6, 1) # 6 variables, order 1
  beam = Beam(d, beta_gamma_ref=beta_gamma_ref, species=species) 
  work = BeamTracking.get_work(beam, Val(n_work))

  track!(beam, ele, work=work)
  p = Particle(beam, 1)

  # 1) Correctness
  @test norm(GTPSA.jacobian(p.v) - M_expected) < tol 
  # 2) Type stability
  if type_stable
    @test_opt track!(beam, ele, work=work)
  end
  # 3) No Allocations
  if no_allocs
    @test @ballocated(track!($beam, $ele, work=$work)) == 0 
  end
end

@testset "BeamTracking.jl" begin
  # Linear tests -----------------------  
  beta_gamma_ref = 1.0
  gamma_ref = sr_gamma(beta_gamma_ref)
  species = Species("electron")
  brho_ref = brho(massof(species), beta_gamma_ref, chargeof(species)) 
  tol = 1e-14 

  # Linear drift test:
  L_d = 5.0
  drift = Linear.Drift(L=L_d)
  M_drift_expected = [1.0  L_d  0.0  0.0  0.0  0.0;
                      0.0  1.0  0.0  0.0  0.0  0.0;
                      0.0  0.0  1.0  L_d  0.0  0.0;
                      0.0  0.0  0.0  1.0  0.0  0.0;
                      0.0  0.0  0.0  0.0  1.0  L_d/gamma_ref^2;
                      0.0  0.0  0.0  0.0  0.0  1.0] 
                      
  test_matrix(drift, 0, M_drift_expected, beta_gamma_ref=beta_gamma_ref, species=species)
  
  # Linear Quadrupole test:
  # Focusing:
  L_q = 1.2
  K1n = 0.36
  B1 = K1n*brho_ref
  qf = Linear.Quadrupole(B1=B1,L=L_q)
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

  test_matrix(qf, 1, M_qf_expected, beta_gamma_ref=beta_gamma_ref, species=species)

  # Defocusing:
  qd = Linear.Quadrupole(B1=-B1,L=L_q)
  M_qd_expected = zeros(6,6)
  M_qd_expected[1:2,1:2] = M_qf_y
  M_qd_expected[3:4,3:4] = M_qf_x
  M_qd_expected[5:6,5:6] = M_qf_z 

  test_matrix(qd, 1, M_qd_expected, beta_gamma_ref=beta_gamma_ref, species=species)
end

#include("linear.jl")
