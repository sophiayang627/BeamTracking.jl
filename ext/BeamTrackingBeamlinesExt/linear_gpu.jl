
function track!(bunch::Bunch, bitsl::AbstractVector{BitsLineElement{S,T,A,B,C}}; work=nothing) where {S,T,A,B,C}
  if length(bitsl) == 0
    return bunch
  end

  # Preallocate the temporaries if not provided
  # This will be slow, the constant dictionary accesses are very type unstable
  # For fast tracking you should preallocate the work array
  if isnothing(work)
    n_temps = 0
    for ele in bitsl
      cur_n_temps = MAX_TEMPS(BeamTracking.ID_TO_TM[ele.tracking_method]) 
      n_temps = cur_n_temps > n_temps ? cur_n_temps : n_temps
    end
    N_particle = get_N_particle(bunch) 
    work = zeros(eltype(bunch.v), N_particle, n_temps)
  end

  # now we have the ENTIRE LATTICE as isbits, so we launch! on this
  launch!(track_isbits!, soaview(bunch), work, bunch.Brho_0, calc_gamma(bunch.species, bunch.Brho_0), bitsl)
  return bunch
end

@inline function track_isbits!(i, v, work, Brho_0, gamma_0, bitsl::AbstractVector{BitsLineElement{S,T,A,B,C}}) where {S,T,A,B,C}
  # We are fully type stable here
  N_ele = length(bitsl)
  @inbounds for j in 1:N_ele
    ele = bitsl[j]
    if ele.tracking_method == BeamTracking.LINEAR_ID
        linear_universal!(i, v, work, Brho_0, gamma_0, ele.L, ele.Brho_ref, ele.BMultipoleParams, ele.BendParams, ele.AlignmentParams)
    else
      error("Tracking method not currently supported for isbits tracking.")
    end
  end
  return
end