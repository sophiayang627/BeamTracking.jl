module Symplectic
using ..BeamTracking: Coords, Beam
using ..GTPSA
using ..AcceleratorLattice

export track!

"""
    track!(ele::Drift, statef::Beam, state0::Beam) -> statef

This function performs exact, hence symplectic, tracking of a `Beam` through a drift.

### Arguments
  - `ele`    -- `Drift`-type element
  - `statef` -- Output `Beam`
  - `state0` -- Input `Beam`
""" track!

function track!(ele::Drift, statef::Beam, statei::Beam)
  @assert !(statef === statei) "Aliasing statef === statei not allowed!"
  L = ele.L
  zi = statei.z
  zf = statef.z

  tilde_m = massof(ele.species_ref) / ele.pc_ref
  β0 = ele.pc_ref / ele.E_tot_ref

  begin
    @. ps = sqrt((1.0 + zi[6])^2 - zi[2]^2 - zi[4]^2)
    @. et = sqrt((1.0 + zi[6])^2 + tilde_m^2)
    @. zf[1] = zi[1] + zi[2] * L / ps
    @. zf[2] = zi[2]
    @. zf[3] = zi[3] + zi[4] * L / ps
    @. zf[4] = zi[4]
    @. zf[5] = zi[5] - (1.0 + zi[6]) * (L / ps - L / (β0 * et))
    @. zf[6] = zi[6]
  end
  return statef
end # function track!(::Drift, ::Beam, ::Beam)


"""
track "linear" quadrupole
"""
function trackQL!(ele::Quadrupole, statef::Beam, statei::Beam)
  @assert !(statef === statei) "Aliasing statef === statei not allowed!"
  L = ele.L
  zi = statei.z
  zf = statef.z

  tilde_m = massof(ele.species_ref) / ele.pc_ref
  β0 = ele.pc_ref / ele.E_tot_ref

  begin
    @. ps = sqrt((1.0 + zi[6])^2 - zi[2]^2 - zi[4]^2)
    @. et = sqrt((1.0 + zi[6])^2 + tilde_m^2)
    @. zf[1] = zi[1] + zi[2] * L / ps
    @. zf[2] = zi[2]
    @. zf[3] = zi[3] + zi[4] * L / ps
    @. zf[4] = zi[4]
    @. zf[5] = zi[5] - (1.0 + zi[6]) * (L / ps - L / (β0 * et))
    @. zf[6] = zi[6]
  end
  return statef
end # function track(Drift)


end # module Symplectic

