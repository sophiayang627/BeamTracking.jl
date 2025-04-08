module BeamTrackingBeamlinesExt
using Beamlines, BeamTracking, GTPSA
using BeamTracking: soaview, get_N_particle, MAX_TEMPS, calc_gamma, launch!
import BeamTracking: track!

# Define my own custom tracking method
# This should be in Beamlines itself so default to this?
struct SciBmadStandard end

function track!(bunch::Bunch, bl::Beamline; work=nothing)
  if length(bl.line) == 0
    return bunch
  end

  # Preallocate the temporaries if not provided
  if isnothing(work)
    n_temps = 0
    for ele in bl.line
      cur_n_temps = MAX_TEMPS(ele.tracking_method) 
      n_temps = cur_n_temps > n_temps ? cur_n_temps : n_temps
    end
    N_particle = get_N_particle(bunch) 
    work = zeros(eltype(bunch.v), N_particle, n_temps)
  end

  for ele in bl.line
    track!(bunch, ele; work=work)
  end

  return bunch
end

function track!(bunch::Bunch, ele::LineElement; work=nothing)
  # Dispatch on the tracking method:
  return track!(bunch, ele, ele.tracking_method; work=work)
end

include("linear.jl")


end