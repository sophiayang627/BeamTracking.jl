isactive(p::AbstractParams) = true
isactive(::Nothing) = false
#isactive(p::BitsBMultipoleParams) = 
#isactive

function linear_universal!(
  i, 
  v, 
  work, 
  Brho_0,
  gamma_0,
  L, 
  Brho_ref, 
  bmultipoleparams, 
  bendparams, 
  alignmentparams
) 
  if isactive(bendparams)
    error("bend tracking not implemented yet")
  end

  if !isactive(bmultipoleparams)
    runkernel!(LinearTracking.linear_drift!, i, v, work, L, L/gamma_0^2)
  else
    if length(bmultipoleparams.bdict) > 1 || !haskey(bmultipoleparams.bdict, 2)
      error("Currently only quadrupole tracking is supported")
    end
    bm1 = bmultipoleparams.bdict[2]
    B1 = bm1.strength
    if bm1.normalized
      B1 *= Brho_ref
    end
    if bm1.integrated
      B1 /= L
    end

    K1 = B1/Brho_0
    mx, my = LinearTracking.linear_quad_matrices(K1, L)
    r56 = L/gamma_0^2
    runkernel!(LinearTracking.linear_coast_uncoupled!, i, v, work, mx, my, r56)
  end
  return v
end

include("linear_cpu.jl")
include("linear_gpu.jl")