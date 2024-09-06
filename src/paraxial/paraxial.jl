"""
Example routine to track through a drift using the paraxial approximation 
and including higher order energy effects. This function can be used both 
for SoA tracking (where `beamf` and `beam0` are the `Beam` struct), and 
also AoS tracking (where `beamf` and `beam0` are just ordinary vectors).
"""
function track!(drift::Drift, beamf, beam0)
  @assert !(beamf === beam0) "Aliasing beamf === beam0 not allowed!"
  
  L = drift.L

  @FastGTPSA begin
  @. zf[1] = z0[1]+z0[2]*L/(1.0+z0[6])
  @. zf[2] = z0[2]
  @. zf[3] = z0[3]+z0[4]*L/(1.0+z0[6])
  @. zf[4] = z0[4]
  @. zf[5] = z0[5]-L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  @. zf[6] = z0[6] 
  end
  return zf
end