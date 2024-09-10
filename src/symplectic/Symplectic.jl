module Symplectic
using ..BeamTracking: Coords, Beam, Coord, Particle
using ..GTPSA
using ..AcceleratorLattice

"""
`track!(ele::Drift, statef::T, state0::T) where {T <: Union{Beam,Particle}} -> statef`

This function performs exact, hence symplectic, tracking of either
a `Particle` or a `Beam` through a drift.

### Arguments
  - `ele`    -- `Drift`-type element
  - `statef` -- Output state, may be either a `Particle` or a `Beam`
  - `state0` -- Input state, may be either a `Particle` or a `Beam`
""" track

function track!(ele::Drift, statef::T, state0::T) where {T <: Union{Beam,Particle}}
  @assert !(statef === state0) "Aliasing statef === state0 not allowed!"
  L = ele.L
  z0 = state0.z
  zf = statef.z

  begin
    @. ps = sqrt((1.0 + z0[6])^2 - z0[2]^2 - z0[4]^2)
    @. et = sqrt((1.0 + z0[6])^2 + tilde_m^2)
    @. zf[1] = z0[1] + z0[2] * L / ps
    @. zf[2] = z0[2]
    @. zf[3] = z0[3] + z0[4] * L / ps
    @. zf[4] = z0[4]
    @. zf[5] = z0[5] - (1.0 + z0[6]) * (L / ps - L / (β0 * et))
    @. zf[6] = z0[6]
  end
  return statef
end # function track(Drift)


end
