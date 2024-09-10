module Paraxial
using ..BeamTracking: Coords, Beam, Coord, Particle
using ..GTPSA
using ..AcceleratorLattice

"""
    track!(ele::Drift, statef::T, state0::T) where {T <: Union{Beam,Particle}} -> statef

Routine to tracking through a drift using the paraial approximation and 
including higher-order energy effects. This function can be used both for 
a single `Particle`, or for a `Beam`.

For SoA tracking, the `Beam` struct should be used, and for `AoS` tracking, 
where `statef` and `state0` are `Vector{Particle{T}}` one can do 
`track!.(ele, statef, state0)` using the dot (`.`) to broadcast.

### Arguments
- `ele`    -- `Drift` type element
- `statef` -- Output state after tracking through, may be either a `Particle` or a `Beam`
- `state0` -- Input state before tracking through, may be either a `Particle` or a `Beam`
"""
function track!(ele::Drift, statef::T, state0::T) where {T <: Union{Beam,Particle}}
  @assert !(statef === state0) "Aliasing statef === state0 not allowed!"
  L = ele.L
  z0 = state0.z
  zf = statef.z

  @FastGTPSA begin
  @. zf[1] = z0[1]+z0[2]*L/(1.0+z0[6])
  @. zf[2] = z0[2]
  @. zf[3] = z0[3]+z0[4]*L/(1.0+z0[6])
  @. zf[4] = z0[4]
  @. zf[5] = z0[5]-L*((z0[2]^2)+(z0[4]^2))/(1.0+z0[6])^2/2.0
  @. zf[6] = z0[6] 
  end
  return statef
end

function helloworld()
  println("hellow world")
end
end

