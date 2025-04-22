
# For BitsBeamline tracking, add the linear tracking method per Beamlines instructions:
function __init__()
  Beamlines.TRACKING_METHOD_MAP[Linear] = 0x1
end
Beamlines.get_tracking_method_extras(::Linear) = SA[]


# Step 1: Unpack the element ---------------------------------------------
function _track!(
  i,
  v,
  work,
  bunch::Bunch,
  ele::Union{LineElement,BitsLineElement}, 
  ::Linear
)
  # Unpack the line element
  ma = ele.AlignmentParams
  bm = ele.BMultipoleParams
  bp = ele.BendParams
  L = ele.L

  # Function barrier
  linear_universal!(i, v, work, bunch, L, bm, bp, ma)
end


# Step 2: Push particles through -----------------------------------------
function linear_universal!(
  i, 
  v, 
  work,
  bunch,
  L, 
  bmultipoleparams, 
  bendparams, 
  alignmentparams
) 
  if isactive(bendparams)
    error("bend tracking not implemented yet")
  end

  gamma_0 = calc_gamma(bunch.species, bunch.Brho_ref)

  if !isactive(bmultipoleparams)
    runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2)
  else 
    if any(t -> t > 2, keys(bmultipoleparams.bdict)) || !haskey(bmultipoleparams.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end

    bm1 = bmultipoleparams.bdict[2]
    K1 = bm1.strength
    if bm1.integrated
      K1 /= L
    end
    if !bm1.normalized
      K1 /= bunch.Brho_ref
    end

    mx, my = LinearTracking.linear_quad_matrices(K1, L)
    r56 = L/gamma_0^2
    runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56)
  end
  return v
end
