module Paraxial
using ..BeamTracking: Coords, Beam
using ..GTPSA
#using ..AcceleratorLattice

struct ParaxialDrift{T}
  L::T
end


"""
    track!(ele::ParaxialDrift, beamf::Beam, beam0::Beam) -> beamf

Routine to tracking through a drift using the paraxial approximation and 
including higher-order energy effects. 

### Arguments
- `ele`   -- `Drift` type element
- `beamf` -- Output beam after tracking through
- `beam0` -- Input beam before tracking through
"""
function track!(L, beamf::Beam, beam0::Beam)
  @assert !(beamf === beam0) "Aliasing beamf === beam0 not allowed!"
  z0 = beam0.z
  zf = beamf.z
  m = beam0.species.mass
  β0 = ele.pc_ref/ ele.E_tot_ref
  
  @FastGTPSA! begin
  @. tilde_m = m*c_light^2/ele.pc_ref*c_light
  @. et = sqrt(1.0 + tilde_m^2/(1+z0[6])^2)
  @. zf[1] = z0[1]+z0[2]*L/(1.0+z0[6])
  @. zf[2] = z0[2]
  @. zf[3] = z0[3]+z0[4]*L/(1.0+z0[6])
  @. zf[4] = z0[4]
  @. zf[5] = z0[5]-L*(1+((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0-1/(β0*et))
  @. zf[6] = z0[6] 
  end
  
  return beamf
end


end
