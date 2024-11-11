module Linear
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

Base.@kwdef struct Drift{T}
  L::T
end

Base.@kwdef struct Quadrupole{T}
  B1::T 
  L::T  
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
    track!(beam::Beam, ele::Linear.Drift) -> beam

Track through a linear drift.

### Arguments
- `beam` -- Input/output beam before/after tracking through
- `ele`  -- `Linear.Drift` type element
"""
function track!(beam::Beam, ele::Linear.Drift; work=get_work(beam, Val{0}()))
  L = ele.L
  v = beam.v
  
  gamma_ref = sr_gamma(beam.beta_gamma_ref)

  @FastGTPSA! begin
  @. v.x  = v.x + L*v.px
  @. v.y  = v.y + L*v.py
  @. v.z  = v.z + L/gamma_ref^2*v.pz
  end

  # Spin unchanged
  
  return beam
end


function track!(beam::Beam, ele::Linear.Quadrupole; work=get_work(beam, Val{1}()))
  v = beam.v
  L = ele.L

  K1n = ele.B1 / brho(massof(beam.species), beam.beta_gamma_ref, chargeof(beam.species))
  gamma_ref = sr_gamma(beam.beta_gamma_ref)

  if K1n >= 0
    k = sqrt(K1n)
    kL = k*L
    sx  = sin(kL)
    cx  = cos(kL)
    sxc = sincu(kL)
    sy  = sinh(kL)
    cy  = cosh(kL)
    syc = sinhcu(kL)
    sgn = 1
  else
    k = sqrt(-K1n)
    kL = k*L
    sx  = sinh(kL)
    cx  = cosh(kL)
    sxc = sinhcu(kL)
    sy  = sin(kL)
    cy  = cos(kL)
    syc = sincu(kL)
    sgn = -1
  end

  
  # Note: 0+work[1] is workaround for GTPSA @FastGTPSA! bug
  @FastGTPSA! begin
  @. work[1]  = cx*v.x + sxc*L*v.px
  @. v.px     = -sgn*k*sx*v.x + cx*v.px
  @. v.x      = 0+work[1] 

  @. work[1]  = cy*v.y + syc*L*v.py
  @. v.py     = sgn*k*sy*v.y + cy*v.py
  @. v.y      = 0+work[1] 

  @. v.z      = v.z + L/gamma_ref^2*v.pz
  end 
 

  return beam
end 

"""
Routine to linearly tracking through a Sector Magnet
"""
function track!(beamf::Beam, ele::Linear.SBend, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.vec
  zf = beamf.vec
  L = ele.L
  q = chargeof(beami.species)

  #curvature of B field
  k = ele.B0 / brho(massof(beami.species),beami.beta_gamma_ref,q)
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
function track!(beamf::Beam, ele::Linear.Combined, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  zi = beami.vec
  zf = beamf.vec
  L = ele.L
  
  brho = brho(massof(beami.species),beami.beta_gamma_ref,q)
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