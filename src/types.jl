#=

Beam, Coord, and Particle type definitions. StructArrays is 
used to handle an internal SoA layout of memory which also 
allows us to mutate, so both the phase space Coord struct 
(defined here) and Quaternion struct (defined in 
ReferenceFrameRotations) can be static. 

Particle struct simply goes from StructArrays SoA to AoS.

=#


# Static phase space coordinate vector
Base.@kwdef struct Coord{T} <: FieldVector{6, T} 
  x::T  = 0.0
  px::T = 0.0
  y::T  = 0.0
  py::T = 0.0
  z::T  = 0.0
  pz::T = 0.0
end

# Static quaternion type defined by ReferenceFrameRotations

struct Beam{T<:StructVector{<:Coord}, U<:Union{Nothing, StructVector{<:Quaternion}}}
  species::Species
  beta_gamma_ref::Float64
  v::T
  q::U
end

# Initialize a Beam with a single particle
function Beam(; species::Species=Species("electron"), beta_gamma_ref=1.0, spin::Union{Bool,Nothing}=nothing,
                x::Number=0.0, px::Number=0.0, y::Number=0.0, py::Number=0.0, z::Number=0.0, pz::Number=0.0)

  T = promote_type(typeof(x), typeof(px), typeof(y),typeof(py), typeof(z), typeof(pz))             
  v = StructArray{Coord{T}}((T[x], T[px], T[y], T[py], T[z], T[pz]))

  if !isnothing(spin)
    if spin == true
      # Use "one(first(..))" etc for GTPSA - descriptor is implicit
      q0 = T[one(first(v.x))]  
      q1 = T[zero(first(v.x))]
      q2 = T[zero(first(v.x))]
      q3 = T[zero(first(v.x))]
      q = StructArray{Quaternion{T}}((q0, q1, q2, q3))
    else
      error("For no spin tracking, please omit the spin kwarg or set spin=nothing. This is to ensure type stability.")
    end
  else
    q = nothing
  end

  return Beam(species, Float64(beta_gamma_ref), v, q)
end


# Initialize Beam as identity GTPSA
# Creates a Beam as identity GTPSA
function Beam(d::Descriptor; species::Species=Species("electron"), beta_gamma_ref=1.0, spin::Union{Bool,Nothing}=nothing)
  GTPSA.numvars(d) == 6 || error("Invalid GTPSA Descriptor! Number of variables must be equal to 6.")
  v = vars(d)
  return Beam(species=species, beta_gamma_ref=beta_gamma_ref, spin=spin, x=v[1], px=v[2], y=v[3], py=v[4], z=v[5], pz=v[6])
end


# Initialize a Beam given some distributions in each coordinate
function Beam(n::Integer; species::Species=Species("electron"), 
              beta_gamma_ref=1.0, spin::Union{Bool,Nothing}=nothing,
              d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
              d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
              d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0) )

  x  = rand(d_x , n)
  px = rand(d_px, n)
  y  = rand(d_y , n)
  py = rand(d_py, n)
  z  = rand(d_z , n)
  pz = rand(d_pz, n)

  v = StructArray{Coord{eltype(x)}}((x, px, y, py, z, pz))
  
  if !isnothing(spin)
    if spin == true
      q0 = ones(eltype(x), n)
      q1 = zeros(eltype(x), n)
      q2 = zeros(eltype(x), n)
      q3 = zeros(eltype(x), n)
      q = StructArray{Quaternion{eltype(x)}}((q0, q1, q2, q3))
    else
      error("For no spin tracking, please omit the spin kwarg or set spin=nothing. This is to ensure type stability.")
    end
  else
    q = nothing
  end

  return Beam(species, Float64(beta_gamma_ref), v, q)
end

struct Particle{T,U<:Union{Nothing,Quaternion{T}}}
  species::Species
  beta_gamma_ref::Float64
  v::Coord{T}
  q::U
end

function Particle(beam::Beam, idx::Integer=1)
  v = beam.v[idx] # StructArrays handles this!
  q = isnothing(beam.q) ? nothing : beam.q[idx]
  return Particle(beam.species, beam.beta_gamma_ref, v, q)
end