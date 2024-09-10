module Paraxial
using ..BeamTracking: Coords, Beam, Coord, Particle
using ..GTPSA
using ..AcceleratorLattice

"""
    track!(L, beamf::Beam, beam0::Beam) -> beamf

Routine to tracking through a drift using the paraxial approximation and 
including higher-order energy effects. 

### Arguments
- `ele`   -- `Drift` type element
- `beamf` -- Output beam after tracking through
- `beam0` -- Input beam before tracking through
"""
function track!(ele::Drift, beamf::Beam, beam0::Beam)
  @assert !(beamf === beam0) "Aliasing beamf === beam0 not allowed!"
  L = ele.L
  z0 = beam0.z
  zf = beamf.z

  @FastGTPSA begin
  @. zf[1] = z0[1]+z0[2]*L/(1.0+z0[6])
  @. zf[2] = z0[2]
  @. zf[3] = z0[3]+z0[4]*L/(1.0+z0[6])
  @. zf[4] = z0[4]
  @. zf[5] = z0[5]-L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  @. zf[6] = z0[6] 
  end
  return beamf
end

end