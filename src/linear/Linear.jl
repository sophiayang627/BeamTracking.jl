module Linear
using ..BeamTracking: Coords, Beam
using ..GTPSA: @FastGTPSA!

export track!

struct Drift{T}
  L::T
end

struct Quadrupole{T}
  L::T
  B1::T
end

struct SBend{T}
  L::T    # Arc length
  B0::T   # Field strength
  g::T    # Coordinate system curvature through element
  e1::T   # Edge 1
  e2::T   # Edge 2
end

struct Solenoid{T}
  L::T
  Bs::T
end

"""
    track!(ele::Linear.Drift, beamf::Beam, beami::Beam) -> beamf

Routine to tracking through a drift using the  approximation and 
including higher-order energy effects. 

### Arguments
- `ele`   -- `Drift` type element
- `beamf` -- Output beam after tracking through
- `beami` -- Input beam before tracking through
"""
function track!(ele::Linear.Drift, beamf::Beam, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  L = ele.L
  zi = beami.z
  zf = beamf.z
  
  @FastGTPSA! begin
  @. zf.x  = zi.x+zi.px*L
  @. zf.px = zi.px
  @. zf.y  = zi.y+zi.py*L
  @. zf.py = zi.py
  @. zf.z  = zi.z-zi.pz*L
  @. zf.pz = zi.pz 
  end
  
  return beamf
end

function track!(ele::Linear.Quadrupole,beamf::Beam,beami::Beam,s::Float64)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.z
  zf = beamf.z
 
  q = beami.species.charge
  pc_ref = sr_pc(q, beami.beta_gamma_0) #reference momentum 
  r = q/pc_ref #rigidity
  k1 = ele.B1 / r #quadrupole strengh

  k = sqrt(abs(k1))
  ks = k * s
  cx = (cosh(ks) * (k1 >= 0) + cos(ks) * (k1 < 0))
  sx = (sinh(ks) * (k1 >= 0) + sin(ks) * (k1 < 0))
  cy = (cos(ks) * (k1 >= 0) + cosh(ks) * (k1 < 0))
  sy = (sin(ks) * (k1 >= 0) + sinh(ks) * (k1 < 0))

  @FastGTPSA! begin
   @. zf.x  = zi.x * cx + zi.px * sx / k 
   @. zf.px = ( -1 * (k1<0) + 1 * (k1>=0) ) zi.x * sx + zi.px * cx
   @. zf.y = zi.y * cy + zi.py * sy
   @. zf.py = ( 1 * (k1<0) - 1 * (k1>=0) )zi.y * k * sy + zi.py * cy
   @. zf.z = zi.z + zi.pz * L 
   @. zf.pz = zi.pz

  end 
end 

end