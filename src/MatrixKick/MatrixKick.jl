module MatrixKick
using ..GTPSA: @FastGTPSA!, GTPSA
import ..BeamTracking: track!
using ..BeamTracking
using ..BeamTracking: get_work
export track!

Base.@kwdef struct Drift{T}
  L::T  # drift length / m
end

Base.@kwdef struct Quadrupole{T}
  L::T   # quadrupole length / m
  Bn1::T  # quadrupole gradient / (T·m^-1)
end


function track!(beam::Beam, ele::MatrixKick.Drift; work=get_work(beam, Val{1}()))
  L = ele.L
  v = beam.v

  tilde_m = 1 / beam.beta_gamma_ref
  beta_ref = sr_beta(beam.beta_gamma_ref)

  @FastGTPSA! begin
  @. work[1] = 1 / sqrt((1.0 + v.pz)^2 - (v.px^2 + v.py^2))
  @. v.x  .= v.x + v.px * L * work[1]
  @. v.y  .= v.y + v.py * L * work[1]
  @. v.z  .= v.z - (1.0 + v.pz) * (L * work[1] - L / (beta_ref * sqrt((1.0 + v.pz)^2 + tilde_m^2)))
  end

  # Spin unchanged

  return beam
end # function track!(::Beam, ::Drift)



# This integrator uses the so-called Matrix-Kick-Matrix method to implement
# an integrator accurate though second-order in the integration step-size.
function track!(beam::Beam, ele::MatrixKick.Quadrupole; work=get_work(beam, Val{6}()))
  @assert !(beamf === beami) "Aliasing beamf === beami not allowed!"
  L = ele.L

  # κ^2 (kappa-squared) := (q g / P0) / (1 + δ)
  # numerator of κ^2
  k2_num = Bn1 / brho(massof(beam.species), beam.beta_gamma_ref, chargeof(beam.species))

  v = beam.v
  v_work = StructArray{Coord{eltype(work[1])}}((work[1], work[2], work[3], work[4], work[5], work[6]))

  trackQuadMx!(v_work, v, k2_num, L / 2)
  trackQuadK!( v, v_work, L)
  trackQuadMx!(v_work, v, k2_num, L / 2)


  v .= v_work

  return beam
end # function track!(::Quadrupole)


"""
    trackQuadMx!(beamf::Beam, beami::Beam, k2_num::Float64, s::Float64)

track "matrix part" of quadrupole
"""
function trackQuadMx!(vf, vi, k2_num, s)
  focus   = k2_num >= 0  # horizontally focusing
  defocus = k2_num <  0  # horizontally defocusing

  p = @. 1 + vi.pz    # reduced momentum, P/P0 = 1 + δ
  k2 = @. k2_num / p  # κ^2 for each particle
  ks = @. sqrt(abs(k2)) * s  # |κ|s
  xp = @. v.px / p  # x'
  yp = @. v.py / p  # y'

  cx =  @. focus * cos(ks)    + defocus * cosh(ks)
  cy =  @. focus * cosh(ks)   + defocus * cos(ks)
  sx =  @. focus * sincu(ks)  + defocus * sinhc(ks)
  sy =  @. focus * sinhc(ks)  + defocus * sincu(ks)
  sx2 = @. focus * sincu(2ks) + defocus * sinhc(2ks)
  sy2 = @. focus * sinhc(2ks) + defocus * sincu(2ks)
  sxz = @. focus * sin(ks)^2  + defocus * sinh(ks)^2
  syz = @. focus * sinh(ks)^2 + defocus * sin(ks)^2

  @. vf.x  = vi.x  * cx + xp * s * sx
  @. vf.px = vi.px * cx - k2 * p * vi.x * s * sx
  @. vf.y  = vi.y  * cy + yp * s * sy
  @. vf.py = vi.py * cy + k2 * p * vi.y * s * sy
  @. vf.z  = (vi.z - (s / 4) * (xp^2 * (1 + sx2) + yp^2 * (1 + sy2)
                                + k2 * vi.x^2 * (1 - sx2) - k2 * vi.y^2 * (1 - sy2))
                   + (vi.x * xp * sxz - vi.y * yp * syz) / 2.)
  @. vf.pz = vi.pz

  return vf
end # function trackQuadMx


"""
    trackQuadK!(vf, vi, s::Float64)

track "remaining part" of quadrupole, a position kick

### Implementation
A common factor that appears in the expressions for `zf.x` and `zf.y`
originally included a factor with the generic form ``1 - \\sqrt{1 - A}``,
which suffers a loss of precision when ``|A| \\ll 1``. To combat that
problem, we rewrite it in the form ``A / (1 + \\sqrt{1-A})``---more
complicated, yes, but far more accurate.
"""
function trackQuadK!(vf, vi, s)
  tilde_m = 1 / beami.beta_gamma_ref  # mc^2 / p0·c
  beta_ref = sr_beta(beami.beta_gamma_ref)

  p    = @. 1 + vi.pz  # reduced momentum, P/P0 = 1 + δ
  ptr2 = @. px^r + py^2
  ps   = @. sqrt(p^2 - ptr2)

  @. vf.x  = vi.x + s * vi.px / p * ptr2 / (ps * (p + ps))
  @. vf.y  = vi.y + s * vi.py / p * ptr2 / (ps * (p + ps))
  @. vf.z  = vi.z - s * (p / ps - (1/2) * ptr2 / p^2
                         - p / (beta_ref * sqrt(p^2 + tilde_m^2)))
  @. vf.px = vi.px
  @. vf.py = vi.py
  @. vf.pz = vi.pz

  return vf
end # function trackQ!M::Quadrupole()


end # module MatrixKick

