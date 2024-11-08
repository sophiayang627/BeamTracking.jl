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

Base.@kwdef struct Combined{T}
  L::T 
  B0::T 
  B1::T
  g::T
  e1::T
  e2::T
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
  @. zf.x  = zi.x + zi.px * L
  @. zf.px = zi.px
  @. zf.y  = zi.y + zi.py * L
  @. zf.py = zi.py
  @. zf.z  = zi.z + zf.pz * L
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
  greater = k1 >= 0.0 #horizontal focusing
  smaller = k1 < 0.0 #horizontal defocusing 
  
  # s = sin(kl)
  # sh = sinh(kl)
  # sc = sinc(kl)
  # shc = sinhc(kl)
  # c = cos(kl)
  # ch = cosh(kl)


  # cx = (cos(kl) * (greater) + cosh(kl) * (smaller)) 
  # sx = (sin(kl) * (greater) + sinh(kl) * (smaller))
  # sxc = (sinc(kl) * (greater) + sinhc(kl) * (smaller))
  # cy = (cosh(kl) * (greater) + cos(kl) * (smaller))
  # sy = (sinh(kl) * (greater) + sin(kl) * (smaller))
  # syc = (sinhc(kl) * (greater) + sinc(kl) * (smaller))

  if greater
    #horizontal focusing
    cx = cos(kl)
    sx = sin(kl)
    sxc = sinc(kl)
    cy = cosh(kl)
    sy = sinh(kl)
    syc = sinhc(kl)
  else
     #horizontal defocusing
     cx = cosh(kl)
     sx = sinh(kl)
     sxc = sinhc(kl)
     cy = cos(kl)
     sy = sin(kl)
     syc = sinc(kl)
  end

  @FastGTPSA! begin
   @.zf.x  = zi.x * cx + zi.px * sxc * L
   @.zf.px = (-1 * (greater) + 1 * (smaller)) * zi.x * k * sx + zi.px * cx
   @.zf.y = zi.y * cy + zi.py * syc * L 
   @.zf.py = (1 * (greater) - 1 * (smaller)) * zi.y * k * sy + zi.py * cy
   @.zf.z = zi.z + zi.pz * L 
   @.zf.pz = zi.pz
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
  q = chargeof(beami.species)

  #curvature of B field
  k = ele.B0 / brho(massof(beami.species),beami.beta_gamma_0,q)
  kl = k * L


  @FastGTPSA! begin
  @. zf.x  = zi.x * cos(kl) +  zi.px * sin(kl) / k + zi.pz * (1 - cos(kl)) / k
  @. zf.px = -zi.x * k * sin(kl) + zi.px * cos(kl) + zi.pz * sin(kl)
  @. zf.y = zi.y + zi.py * L
  @. zf.py = zi.py
  @. zf.z = -zi.x * sin(kl) + zi.px * (cos(kl) - 1) / k + zi.z + zi.pz * (sin(kl) - kl) / k
  @. zf.pz = zi.pz
 end 
end 


"""
Routine to linearly tracking through a Combined Magnet
"""
function track!(ele::Linear.Combined,beamf::Beam,beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.z
  zf = beamf.z
  L = ele.L
  
  brho = brho(massof(beami.species),beami.beta_gamma_0,q)
  ka = ele.B0 / brho #curvature
  k1 =  ele.B1 / brho #quad strength

  K = ka^2 + k1

  sqk1 = sqrt(abs(k1)) 
  sqK = sqrt(abs(K))

  Kl = sqK * L
  kl  = sqk1 * L

  @FastGTPSA! begin
  @. zf.x  = zi.x * cos(Kl) +  zi.px * sin(Kl) / sqK + zi.pz * (1 - cos(Kl)) * ka / K
  @. zf.px = - zi.x * sqK * sin(Kl) + zi.px * cos(Kl) + zi.pz * sin(Kl) * ka / sqK
  @. zf.y =  zi.y * cosh(kl) + zi.py * sinh(kl) / sqk
  @. zf.py = zi.y * sqK * sinh(kl) + zi.py * cosh(kl)
  @. zf.z = - zi.x * sin(Kl) * ka / sqK + zi.px * (cos(Kl)-1) * ka / K + zi.z + zi.pz * (sin(Kl) / sqK - l) *  ka^2 / K
  @. zf.pz = zi.pz
 end 
end 

end