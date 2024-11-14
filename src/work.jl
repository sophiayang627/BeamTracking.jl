"""
    get_work(bunch::Bunch, ::Val{N}) where N -> work

Returns a tuple of `N` arrays of type `eltype(Bunch.v.x)` and 
length `length(Bunch.v.x)` which may be used as temporaries.

### Arguments
- `bunch`     -- Bunch to extract types and number of particles from
- `::Val{N}` -- Number of `N` temporary arrays desired
"""
function get_work(bunch::Bunch, ::Val{N}) where {N}
  sample = first(bunch.v.x)
  T = typeof(sample)
  N_particle = length(bunch.v.x)

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