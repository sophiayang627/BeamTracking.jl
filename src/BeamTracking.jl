module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays,
      Distributions,
      Unitful


include("aapc.jl")

export Beam,
       Coords,
       Particle,
       Coord, 
       MatrixKick,
       Linear,
       sr_gamma, 
       sr_gamma_m1,
       sr_beta,
       sr_pc,
       sr_ekin,
       sr_etot,
       brho,
       chargeof,
       massof,
       Species,
       sincu,
       sinhc,
       track!

# SoA ----------------------------------
Base.@kwdef struct Coords{T}
  x::Vector{T}
  px::Vector{T}
  y::Vector{T}
  py::Vector{T}
  z::Vector{T}
  pz::Vector{T}
end

function Coords(
  n::Integer;
  d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
  d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
  d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0)
)
  x  = rand(d_x , n)
  px = rand(d_px, n)
  y  = rand(d_y , n)
  py = rand(d_py, n)
  z  = rand(d_z , n)
  pz = rand(d_pz, n)

  return Coords(x, px, y, py, z, pz)
end



struct Beam{S,T}
  species::Species
  beta_gamma_0::S
  z::Coords{T}
end

function Beam(
  n::Integer; species::Species=Species("electron"), beta_gamma_0=1,
  d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
  d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
  d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0),
)

  coords = Coords(n; d_x=d_x, d_px=d_px, d_y=d_y, d_py=d_py, d_z=d_z, d_pz=d_pz)

  return Beam(species, beta_gamma_0, coords)
end

# Creates a Beam as identity GTPSA
function Beam(d::Descriptor; species::Species=Species("electron"), beta_gamma_0=1)
  GTPSA.numvars(d) == 6 || error("Invalid GTPSA Descriptor! Number of variables must be equal to 6.")
  z = vars(d)
  return Beam(species, beta_gamma_0, Coords([z[1]], [z[2]], [z[3]], [z[4]], [z[5]], [z[6]]))
end

#Create a Beam with single particle for testing
function Beam(
  ;species::Species=Species("electron"), beta_gamma_0=1,
  x=0.0, px=0.0, y=0.0, py=0.0, z=0.0, pz=0.0,
  )
  
  coords = Coords(x=[x], px=[px], y=[y], py=[py], z=[z], pz=[pz])

  return Beam(species, beta_gamma_0, coords)
end

#AoS from SoA for single particle
# Extract the phase space coord of a particle in a beam 
Base.@kwdef struct Coord{T} <: FieldVector{6, T} # Just like Coords but only 1 Coord 
  x::T  = 0.0
  px::T = 0.0
  y::T  = 0.0
  py::T = 0.0
  z::T  = 0.0
  pz::T = 0.0
end

struct Particle{S,T}
  species::Species
  beta_gamma_0::S
  z::Coord{T}
end

function Particle(b::Beam, n::Integer=1)
  z = b.z
  coord = Coord(z.x[n],z.px[n],z.y[n],z.py[n],z.z[n],z.pz[n])
  
  return Particle(b.species, b.beta_gamma_0, coord)
end

# Empty tracking method ----------------
track!(beamf::Beam, ::Nothing, beami::Beam) = (beamf.z .= beami.z; return)

# --------------------------------------
include("utils.jl")

# Modules separated:
include("MatrixKick/MatrixKick.jl") 
include("Linear/Linear.jl")   


end
