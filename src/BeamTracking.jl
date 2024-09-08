module BeamTracking
using AcceleratorLattice, 
      GTPSA,
      ReferenceFrameRotations,
      StaticArrays

# Temporary until AtomicAndPhysicalConstants is cleaned up ----
import AtomicAndPhysicalConstants: AtomicAndPhysicalConstants
const Species = AtomicAndPhysicalConstants.Particle
# -------------------------------------------------------------

export Beam, 
       Symplectic,
       Paraxial

# SoA ----------------------------------
struct Coords{T} <: FieldVector{6, T}
  x::Vector{T}
  px::Vector{T}
  y::Vector{T}
  py::Vector{T}
  z::Vector{T}
  pz::Vector{T}
end

struct Beam{T} # Must agree exactly with `Particle`!
  species::Species
  z::Coords{T}
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

struct Particle{T} # Must agree exactly with `Beam`!
  species::Species
  z::Coord{T}      
end
# ---------------------------------------


# Modules separated:
include("symplectic/Symplectic.jl") 
include("paraxial/Paraxial.jl")    

end