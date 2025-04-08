
function track!(bunch::Bunch, bl::Beamline; work=nothing, static=true)
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

  if !static
    for ele in bl.line
      track!(bunch, ele; work=work)
    end
  else # get the entire lattice as one track_chain
    chain = make_track_chain(bunch, bl)
    launch!(chain, bunch.v, work)
  end

  return bunch
end



function make_track_chain(bunch, bl::Beamline)
  chain = KernelCall[]
  for ele in bl.line
    append!(chain, make_track_chain(bunch, ele, ele.tracking_method))
  end
  return chain
end

#=
function track!(bunch::Bunch, ele::LineElement; work=zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(ele.tracking_method)))
  # Build a track chain for this tracking method:
  track_chain = make_track_chain(bunch, ele, ele.tracking_method)[1]
  
  # Now launch the tracking
  launch!(track_chain[1], bunch.v, work, track_chain[2:end]...)
  return bunch
end
=#
#=function track!(
  bunch::Bunch, 
  ele::LineElement, 
  ::Linear; 
  work=zeros(eltype(bunch.v), get_N_particle(bunch), MAX_TEMPS(Linear())),
)#=
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L
  Brho_ref = ele.Brho_ref
=#
  # Function barrier to build the tracking chain for this element
  #track_chain = linear_track_chain(bunch, ma, bm, bp, L, Brho_ref)

  track_chain = make_track_chain(bunch, ele, Linear())

  # Now launch
  launch!(track_chain, bunch.v, work)
  return bunch
  
  # Function barrier (this function is now fully compiled)
  # For some reason, inlining this is faster/zero allocs
  # copy-paste is slower and so is @noinline so I guess LLVM is 
  # doing some kind of constant propagation while inlining this?
  #return @inline _track_linear!(bunch, ma, bm, bp, L, Brho_ref, work)
end
=#
function make_track_chain(bunch::Bunch, ele::LineElement, ::Linear)
    # Unpack the line element
    ma = ele.AlignmentParams
    bm = ele.BMultipoleParams
    bp = ele.BendParams
    L = ele.L
    Brho_ref = ele.Brho_ref

    # Function barrier to build the tracking chain for this element
    track_chain = linear_track_chain(bunch, ma, bm, bp, L, Brho_ref)
    return track_chain
end


@inline function linear_track_chain(
  bunch::Bunch, 
  ma::Union{AlignmentParams,Nothing},
  bm::Union{BMultipoleParams,Nothing}, 
  bp::Union{BendParams,Nothing},
  L, 
  Brho_ref
) 

  chain = KernelCall[]
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_0)

  if !isnothing(bp)
    error("Bend tracking not implemented yet")
  end

  if !isnothing(ma)
    push!(chain, KernelCall(ExactTracking.misalign!, (ma.x_offset, ma.y_offset, -1)))
  end

  if isnothing(bm) || length(bm.bdict) == 0 # Drift
    #launch!(LinearTracking.linear_drift!, v, nothing, L, L/gamma_0^2)
    push!(chain, KernelCall(LinearTracking.linear_drift!, (L, L/gamma_0^2)))
  else
    if length(bm.bdict) > 1 || !haskey(bm.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end

    bm1 = bm.bdict[2]
    if bm1.tilt != 0
      error("Currently multipole tilts not supported")
    end

    # Tracking should work with B-fields. The lattice has no understanding of a particle beam
    # or the different types of particles shooting through

    # The Bunch struct itself can store the normalization factor of its own phase space coordinates 
    # but that may in general be different from the lattice.

    # So we always just get the B-fields from the lattice
    B1 = bm1.strength
    if bm1.normalized
      B1 *= Brho_ref
    end
    if bm1.integrated
      B1 /= L
    end

    # Now get K1 from the bunch itself:
    K1 = B1/bunch.Brho_0

    mx, my = LinearTracking.linear_quad_matrices(K1, L)
    r56 = L/gamma_0^2 
    #launch!(LinearTracking.linear_coast_uncoupled!, v, wq, mx, my, r56)
    push!(chain, KernelCall(LinearTracking.linear_coast_uncoupled!, (mx, my, r56)))
  end

  if !isnothing(ma)
    push!(chain, KernelCall(ExactTracking.misalign!, (ma.x_offset, ma.y_offset, -1)))
  end
  
  return chain
end

#=

function _track_linear!(
  bunch::Bunch, 
  ma::Union{AlignmentParams,Nothing},
  bm::Union{BMultipoleParams,Nothing}, 
  bp::Union{BendParams,Nothing},
  L, 
  Brho_ref;
  work=nothing,
) 

  chain = ()
  v = soaview(bunch)
  gamma_0 = calc_gamma(bunch.species, bunch.Brho_0)

  if !isnothing(bp)
    error("Bend tracking not implemented yet")
  end

  if !isnothing(ma)
    chain = (chain, (ExactTracking.misalign!, ma.x_offset, ma.y_offset, -1))
  end

  if isnothing(bm) || length(bm.bdict) == 0 # Drift
    #launch!(LinearTracking.linear_drift!, v, nothing, L, L/gamma_0^2)
    chain = (chain, (LinearTracking.linear_drift!, L, L/gamma_0^2))
  else
    if length(bm.bdict) > 1 || !haskey(bm.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end

    bm1 = bm.bdict[2]
    if bm1.tilt != 0
      error("Currently multipole tilts not supported")
    end

    # Tracking should work with B-fields. The lattice has no understanding of a particle beam
    # or the different types of particles shooting through

    # The Bunch struct itself can store the normalization factor of its own phase space coordinates 
    # but that may in general be different from the lattice.

    # So we always just get the B-fields from the lattice
    B1 = bm1.strength
    if bm1.normalized
      B1 *= Brho_ref
    end
    if bm1.integrated
      B1 /= L
    end

    # Now get K1 from the bunch itself:
    K1 = B1/bunch.Brho_0

    mx, my = LinearTracking.linear_quad_matrices(K1, L)
    r56 = L/gamma_0^2 
    work = isnothing(work) ? zeros(eltype(bunch.v), get_N_particle(bunch), 1) : work
    #launch!(LinearTracking.linear_coast_uncoupled!, v, wq, mx, my, r56)
    chain = (chain, (LinearTracking.linear_coast_uncoupled!, mx, my, r56))
  end

  if !isnothing(ma)
    chain = (chain, (ExactTracking.misalign!, ma.x_offset, ma.y_offset, -1))
  end
  
  launch!(chain[2:end], v, work)
  return bunch
end=#