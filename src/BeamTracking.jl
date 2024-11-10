module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays,
      StructArrays,
      Distributions,
      Unitful


include("aapc.jl")

export Beam,       
       Coord, 
       Quaternion
       Particle,
       
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





include("types.jl")

# Empty tracking method to be imported by submodules 
track!(beam::Beam, ::Nothing) = beam

# --------------------------------------------------

include("utils.jl")

# Modules separated:
include("MatrixKick/MatrixKick.jl") 
include("Linear/Linear.jl")   
#include("Misc/Misc.jl")


end
