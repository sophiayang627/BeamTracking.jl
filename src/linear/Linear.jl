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
  zi = beami.z
  zf = beamf.z
  
  @FastGTPSA! begin
  @. zf.x  = zi.x+zi.px*L
  @. zf.px = zi.px
  @. zf.y  = zi.y+zi.py*L
  @. zf.py = zi.py
  @. zf.z  = zi.z-L*((zi.px^2)+(zi.py^2))/(1.0+zi.pz)^2/2.0
  @. zf.pz = zi.pz 
  end
  
  return beamf
end


end