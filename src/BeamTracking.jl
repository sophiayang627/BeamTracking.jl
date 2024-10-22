module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays,
      Distributions


include("aapc.jl")

export Beam,
       Coords,
       Symplectic,
       Linear,
        sr_gamma, 
        sr_gamma_m1,
        sr_beta,
        sr_pc,
        sr_ekin,
        sr_etot,
        brho

# SoA ----------------------------------
Base.@kwdef struct Coords{T} <: FieldVector{6, T}
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

#Create a Beam struct with single particle for testing
function Beam(;species::Species=Species("electron"), beta_gamma_0=1,
  x::Float64=0.0, px::Float64=0.0, y::Float64=0.0,
  py::Float64=0.0, z::Float64=0.0,pz::Float64=0.0,
  )
  coords = Coords(x=[x], px=[px], y=[y], py=[py], z=[z], pz=[pz])
  return Beam(species, beta_gamma_0, coords)
end

# --------------------------------------
include("utils.jl")

# Modules separated:
include("symplectic/Symplectic.jl") 
include("linear/Linear.jl")    

end
