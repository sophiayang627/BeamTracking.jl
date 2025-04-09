
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

function track!(bunch::Bunch, ele; work=nothing)
  # Dispatch on the tracking method:
  return track!(bunch, ele, ele.tracking_method, work)
end


function track!(
  bunch::Bunch, 
  ele::LineElement, 
  tm::Linear,
  work=zeros(eltype(bunch.v), N_particle, MAX_TEMPS(tm))
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L
  Brho_ref = ele.Brho_ref

  # Function barrier (this function is now fully compiled)
  # For some reason, inlining this is faster/zero allocs
  # copy-paste is slower and so is @noinline so I guess LLVM is 
  # doing some kind of constant propagation while inlining this?
  return _track_linear!(bunch, ma, bm, bp, L, Brho_ref, work)
end

function _track_linear!(bunch, ma, bm, bp, L, Brho_ref, work)
  v = soaview(bunch)
  linear_universal!(nothing, v, work, bunch.Brho_0, calc_gamma(bunch.species, bunch.Brho_0), L, Brho_ref, bm, bp, ma)
  return bunch
end
