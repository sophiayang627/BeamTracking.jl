module Linear
using ..GTPSA: @FastGTPSA!, GTPSA
using ..BeamTracking


export track!

Base.@kwdef struct Drift{T}
  L::T
end

Base.@kwdef struct Quadrupole{T}
  L::T
  B1::T
end

Base.@kwdef struct SBend{T}
  L::T    # Arc length
  B0::T   # Field strength
  g::T    # Coordinate system curvature through element
  e1::T   # Edge 1
  e2::T   # Edge 2
end

Base.@kwdef struct Solenoid{T}
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


"""
Routine to linearly tracking through a quadrupole
"""
function track!(ele::Linear.Quadrupole,beamf::Beam,beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.z
  zf = beamf.z
  L = ele.L
  q = chargeof(beami.species)

  k1 = ele.B1 / brho(massof(beami.species),beami.beta_gamma_0,q) #quadrupole strengh
  k = sqrt(abs(k1))
  kl = k * L
  greater = k1 >= 0 #horizontal focusing
  smaller = k1 < 0 #vertical focusing 

  cx = (cos(kl) * (greater) + cosh(kl) * (smaller)) 
  sx = (sin(kl) * (greater) + sinh(kl) * (smaller))
  cy = (cosh(kl) * (greater) + cos(kl) * (smaller))
  sy = (sinh(kl) * (greater) + sin(kl) * (smaller))

  @FastGTPSA! begin
   @. zf.x  = zi.x * cx + zi.px * sx / k 
   @. zf.px = ( -1 * (smaller) + 1 * (greater) ) * zi.x * sx + zi.px * cx
   @. zf.y = zi.y * cy + zi.py * sy
   @. zf.py = ( 1 * (smaller) - 1 * (greater) ) * zi.y * k * sy + zi.py * cy
   @. zf.z = zi.z + zi.pz * L 
   @. zf.pz = zi.pz
  end 

  return beamf
end 

"""
Routine to linearly tracking through a Sector Magnet
"""
function track!(ele::Linear.SBend,beamf::Beam,beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.z
  zf = beamf.z
  L = ele.L

end 



end