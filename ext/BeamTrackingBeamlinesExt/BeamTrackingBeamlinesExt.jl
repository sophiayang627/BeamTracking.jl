module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA, StaticArrays
using Beamlines: AbstractBitsParams, BitsBMultipoleParams, BitsLineElement, BitsBendParams, BitsAlignmentParams
using BeamTracking: soaview, get_N_particle, MAX_TEMPS, calc_gamma, launch!, runkernel!
import BeamTracking: track!, make_track_chain

function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.tracking_method_extras(::Linear) = SA[]

# Define my own custom tracking method
# This should be in Beamlines itself so default to this?
struct SciBmadStandard end

# We now define the Linear tracking method per instructions in Beamlines:



include("linear.jl")


end