module MatrixKick
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

struct Drift{T}
  L::T
end

struct Quadrupole{T}
  L::T
  B1::T
end

"""
    track!(beamf::Beam, ele::Drift, beami::Beam) -> beamf

This function performs exact, hence symplectic, tracking of a `Beam` through a drift.

### Arguments
  - `beamf` -- Output `Beam`
  - `ele`    -- `Drift`-type element
  - `beami` -- Input `Beam`
"""
function track!(beamf::Beam, ele::Drift, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  @assert !(beamf.species === beami.species) "Input species must be equal to output!"
  L = ele.L
  zi = beami.vec
  zf = beamf.vec

  tilde_m = mass(beami.species) / pc_ref(beami.species, beami.beta_0)

  @. ps = sqrt((1.0 + zi.pz)^2 - zi.px^2 - zi.py^2)
  @. et = sqrt((1.0 + zi.pz)^2 + tilde_m^2)
  @. zf.x = zi.x + zi.px * L / ps
  @. zf.px = zi.px
  @. zf.y = zi.y + zi.py * L / ps
  @. zf.py = zi.py
  @. zf.z = zi.z - (1.0 + zi.pz) * (L / ps - L / (beami.beta_0 * et))
  @. zf.pz = zi.pz

  return beamf
end # function track!(::Drift, ::Beam, ::Beam)


"""
track quadrupole
"""
function track!(beamf::Beam, ele::Quadrupole, beami::Beam)
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  L = ele.L
  zi = beami.vec
  zf = beamf.vec

  tilde_m = massof(ele.species_ref) / ele.pc_ref
  β0 = ele.pc_ref / ele.E_tot_ref

  @. ps = sqrt((1.0 + zi.pz)^2 - zi.px^2 - zi.py^2)
  @. et = sqrt((1.0 + zi.pz)^2 + tilde_m^2)
  @. zf.x = zi.x + zi.px * L / ps
  @. zf.px = zi.px
  @. zf.y = zi.y + zi.py * L / ps
  @. zf.py = zi.py
  @. zf.z = zi.z - (1.0 + zi.pz) * (L / ps - L / (β0 * et))
  @. zf.pz = zi.pz

  return beamf
end # function track!(::Quadrupole)

"""
track "matrix part" of quadrupole
"""
function trackQuadMf!(beamf::Beam, beami::Beam, s::Float64, kappa_num::Float64)
  zi = beami.vec
  zf = beamf.vec

  tilde_m = mass(beami.species) / pc_ref(beami.species, beami.beta_0)

  p_red = 1 + zi.pz # reduced Pz

  κ = kappa_num / p_red
  κs = κ * s
  cx = (cos(κs) * (k1 >= 0) + cosh(ks)*(k1 < 0))
  sx = sin(κs)
  cy = cosh(κs)
  sy = sinh(κs)
  sx2 = sinc(2κs)
  sy2 = sinch(2κs)
  xp = px / p_red
  yp = py / p_red
  
  @. zf.x  = zi.x * cx + xp * sx / κ
  @. zf.px = zi.px * cx - p_red * κ * zi.x
  @. zf.y  = zi.y * cy + yp * sy / κ
  @. zf.py = zi.py * cy + p_red * κ * zi.y
  @. zf.z  = (zi.z - (s / 4.) * (xp^2 * (1. + sx2) + yp^2 * (1 + sy2)
                                  + (κ * x)^2 * (1 - sx2) + (κ * y)^2 * (1 - sy2))
                    + (x * xp * sx^2 - y * yp * sy^2) / 2.)
  @. zf.pz = zi.pz

  return beamf
end # function trackQ!M::Quadrupole()


end # module MatrixKick

