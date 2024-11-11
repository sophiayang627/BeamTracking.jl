"""
    get_work(beam::Beam, ::Val{N}) where N -> work

Returns a tuple of `N` arrays of type `eltype(Beam.v.x)` and 
length `length(Beam.v.x)` which may be used as temporaries.

### Arguments
- `beam`     -- Beam to extract types and number of particles from
- `::Val{N}` -- Number of `N` temporary arrays desired
"""
function get_work(beam::Beam, ::Val{N}) where {N}
  sample = first(beam.v.x)
  T = typeof(sample)
  N_particle = length(beam.v.x)

  # Instead of using zeros, we do this to ensure 
  # same GTPSA descriptor if T isa TPS.
  return ntuple(Val{N}()) do t
    r = Vector{T}(undef, N_particle)
    for idx in eachindex(r)
      r[idx] = zero(sample)
    end
    r
  end
end