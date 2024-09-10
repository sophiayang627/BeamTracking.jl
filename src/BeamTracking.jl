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

struct Beam{T}
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


# Modules separated:
include("symplectic/Symplectic.jl") 
include("paraxial/Paraxial.jl")    

end