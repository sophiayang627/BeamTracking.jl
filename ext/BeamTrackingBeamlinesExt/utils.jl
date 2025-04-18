
function get_work(bunch::Bunch, bl::Beamline)
  n_temps = 0
  for ele in bl.line
    cur_n_temps = MAX_TEMPS(ele.tracking_method) 
    n_temps = cur_n_temps > n_temps ? cur_n_temps : n_temps
  end
  N_particle = get_N_particle(bunch) 
  work = zeros(eltype(bunch.v), N_particle, n_temps)
  return work
end

function get_work(bunch::Bunch, bbl::BitsBeamline{TM}) where {TM}
  # Preallocate the temporaries if not provided
  # This will be slow, the constant dictionary accesses are very type unstable
  # For fast tracking you should preallocate the work array
  if TM == Beamlines.MultipleTrackingMethods
    INVERSE_TRACKING_METHOD_MAP = Dict(value => key for (key, value) in Beamlines.TRACKING_METHOD_MAP)
    n_temps = 0
    for i in 1:length(bbl.tracking_method)
      tm = bbl.tracking_method[i]
      tme = bbl.tracking_method_extras[i]
      return INVERSE_TRACKING_METHOD_MAP
      cur_n_temps = MAX_TEMPS(INVERSE_TRACKING_METHOD_MAP[tm](tme...)) 
      n_temps = cur_n_temps > n_temps ? cur_n_temps : n_temps
    end
  else
    n_temps = MAX_TEMPS(TM)
  end
  N_particle = get_N_particle(bunch) 
  work = zeros(eltype(bunch.v), N_particle, n_temps)
  return work
end


function check_Brho(Brho_ref, bunch::Bunch)
  if !isnan(Brho_ref) && !(Brho_ref ≈ bunch.Brho_0)
    @warn "The reference energy of the bunch does NOT equal the reference energy of the Beamline. 
    Normalized field strengths in tracking always use the reference energy of the bunch. "
  end
end
