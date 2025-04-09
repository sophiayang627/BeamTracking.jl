module BeamTracking
using GTPSA,
      ReferenceFrameRotations,
      StaticArrays, 
      SIMD,
      VectorizationBase,
      Polyester,
      LoopVectorization
      
import GTPSA: sincu, sinhcu
import Base: setproperty!

export Bunch, Species, ParticleView, ELECTRON, POSITRON, PROTON, ANTIPROTON, sincu, sinhcu
export LinearTracking, Linear
export ExactTracking, Exact
export track!, make_track_chain, KernelCall

include("utils.jl")
include("kernel.jl")
include("types.jl")

# Stores the objectid of each tracking method to its type
const ID_TO_TM = Dict{UInt, Any}()

include("modules/ExactTracking.jl") #; TRACKING_METHOD(::ExactTracking) = Exact
include("modules/LinearTracking.jl") #; TRACKING_METHOD(::LinearTracking) = Linear


# Empty tracking method to be imported+implemented by package extensions
function track! end

function make_track_chain end

function MAX_TEMPS end
# --------------------------------------------------


# Modules separated:
#include("MatrixKick/MatrixKick.jl")
#include("Linear/Linear.jl")


end
