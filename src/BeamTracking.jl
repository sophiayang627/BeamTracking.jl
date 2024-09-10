module BeamTracking
using AcceleratorLattice, 
      GTPSA,
      ReferenceFrameRotations,
      StaticArrays,
      Distributions

# Temporary until AtomicAndPhysicalConstants is cleaned up ----
import AtomicAndPhysicalConstants: AtomicAndPhysicalConstants
const Species = AtomicAndPhysicalConstants.Particle
# -------------------------------------------------------------

export Beam,
       Coords,
       Coord,
       Particle,
       Symplectic,
       Paraxial,
       Species

# SoA ----------------------------------
struct Coords{T} <: FieldVector{6, T}
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

struct Beam{T} # Must agree exactly with `Particle`!
  species::Species
  z::Coords{T}
end

function Beam(
  n::Integer; species::Species=Species("electron"),
  d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
  d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
  d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0)
)

  coords = Coords(n; d_x=d_x, d_px=d_px, d_y=d_y, d_py=d_py, d_z=d_z, d_pz=d_pz)

  return Beam(species, coords)
end

# --------------------------------------

# AoS ----------------------------------
struct Coord{T} <: FieldVector{6, T}
  x::T
  px::T
  y::T
  py::T
  z::T
  pz::T
end

function Coord(;
  d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
  d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
  d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0)
)
  x  = rand(d_x )
  px = rand(d_px)
  y  = rand(d_y )
  py = rand(d_py)
  z  = rand(d_z )
  pz = rand(d_pz)

  return Coord(x, px, y, py, z, pz)
end

struct Particle{T} # Must agree exactly with `Beam`!
  species::Species
  z::Coord{T}      
end

function Particle(;
  species::Species=Species("electron"),
  d_x::Distribution=Normal(0,0), d_px::Distribution=Normal(0,0), 
  d_y::Distribution=Normal(0,0), d_py::Distribution=Normal(0,0), 
  d_z::Distribution=Normal(0,0), d_pz::Distribution=Normal(0,0)
)

  coord = Coord(d_x=d_x, d_px=d_px, d_y=d_y, d_py=d_py, d_z=d_z, d_pz=d_pz)

  return Particle(species, coord)
end

# ---------------------------------------


# Modules separated:
include("symplectic/Symplectic.jl") 
include("paraxial/Paraxial.jl")    

end