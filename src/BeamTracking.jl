module BeamTracking
using AcceleratorLattice, 
      GTPSA,
      ReferenceFrameRotations,
      StaticArrays

# Temporary until AtomicAndPhysicalConstants is cleaned up ====
import AtomicAndPhysicalConstants: AtomicAndPhysicalConstants
const Species = AtomicAndPhysicalConstants.Particle
# =============================================================

export Beam, 
       Symplectic,
       Paraxial

struct Beam{T} <: FieldVector{6, T}
  x::Vector{T}
  px::Vector{T}
  y::Vector{T}
  py::Vector{T}
  z::Vector{T}
  pz::Vector{T}
end

# Modules fully separated:

module Symplectic
using ..BeamTracking: Beam
using ..GTPSA
using ..AcceleratorLattice
include("symplectic/symplectic.jl") # If more files added, add more includes here!
end

module Paraxial
using ..BeamTracking: Beam
using ..GTPSA
using ..AcceleratorLattice
include("paraxial/paraxial.jl")     # If more files added, add more includes here!
end

end