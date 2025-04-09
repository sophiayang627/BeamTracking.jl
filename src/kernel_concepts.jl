
const REGISTER_SIZE = VectorizationBase.register_size()
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6



struct KernelCall{F<:Function,T<:Tuple}
  f!::F
  args::T
end


@inline function launch!(
  chain::C,
  v::A,
  work;
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {C<:Vector{<:KernelCall},A}
  N_particle = size(v, 1)
  N_fc = length(chain)
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    if N_particle >= multithread_threshold
      Threads.@threads for j in 1:simd_lane_width:N_SIMD
        i = lane+j
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        for k in 1:N_fc
          @inbounds f! = chain[k].f!
          @inbounds args = chain[k].args
          f!(i, v, work, args...)
        end
      end
    else
      for j in 1:simd_lane_width:N_SIMD
        i = lane+j
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        for k in 1:N_fc
          @inbounds f! = chain[k].f!
          @inbounds args = chain[k].args
          f!(i, v, work, args...)
        end
      end
    end
    # Do the remainder
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!" 
      for k in 1:N_fc
        @inbounds f! = chain[k].f!
        @inbounds args = chain[k].args
        f!(i, v, work, args...)
      end
    end
  else
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        for k in 1:N_fc
          @inbounds f! = chain[k].f!
          @inbounds args = chain[k].args
          f!(i, v, work, args...)
        end
      end
    else
      @simd for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        for k in 1:N_fc
          @inbounds f! = chain[k].f!
          @inbounds args = chain[k].args
          f!(i, v, work, args...)
        end
      end
    end
  end
  return v
end













function recursive_unroll(chain::Tuple, i, v, work)
  if length(chain) == 0
    return
  end
  f! = chain[1][1]
  args = Base.tail(chain[1])
  f!(i, v, work, args...)
  recursive_unroll(Base.tail(chain), i, v, work)
end

@inline function launch!(
  chain::C,
  v::A,
  work;
  simd_lane_width=16, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {C<:NTuple,A}
  N_particle = size(v, 1)
  N_fc = length(chain)
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    if N_particle >= multithread_threshold
      Threads.@threads for j in 1:simd_lane_width:N_SIMD
        i = lane+j
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        for k in 1:N_fc
          @inbounds f! = chain[k][1]
          @inbounds args = chain[k][2:end]
          f!(i, v, work, args...)
        end
      end
    else
      for j in 1:simd_lane_width:N_SIMD
        i = lane+j
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        for k in 1:N_fc
          @inbounds f! = chain[k][1]
          @inbounds args = chain[k][2:end]
          f!(i, v, work, args...)
        end
      end
    end
    # Do the remainder
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!" 
      for k in 1:N_fc
        @inbounds f! = chain[k][1]
        @inbounds args = chain[k][2:end]
        f!(i, v, work, args...)
      end
    end
  else
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        for k in 1:N_fc
          @inbounds f! = chain[k][1]
          @inbounds args = chain[k][2:end]
          f!(i, v, work, args...)
        end
      end
    else
      @simd for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        for k in 1:N_fc
          @inbounds f! = chain[k][1]
          @inbounds args = chain[k][2:end]
          f!(i, v, work, args...)
        end
      end
    end
  end
  return v
end



# Generic function to launch a kernel on the bunch coordinates matrix
# Matrix v should ALWAYS be in SoA whether for real or as a view via tranpose(v)

"""
    launch!(f!::F, v, v0, work, args...; simd_lane_width, multithread_threshold)

General purpose function to launch a kernel `f!`. The syntax for a kernel `f!` must 
ALWAYS be the following:

## Arguments
- `i`       -- Particle index
- `v`       -- Input/output matrix as an SoA or SoA view ALWAYS! (use transpose if AoS)
- `work`    -- A Vector of temporary vectors (columns of v) to run the kernel `f!`
- `args...` -- Any further arguments to run the kernel

## Keyword Arguments
- `simd_lane_width`       -- The number of SIMD lanes to use. Default is `REGISTER_SIZE/sizeof(eltype(A))`
- `multithread_threshold` -- Number of particles at which multithreading is used. Default is `1e6``
"""
@inline function launch!(
  f!::F, 
  v::A,
  work, 
  args...; 
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {F<:Function,A}
  N_particle = size(v, 1)
  if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
    lane = VecRange{simd_lane_width}(0)
    rmn = rem(N_particle, simd_lane_width)
    N_SIMD = N_particle - rmn
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    else
      for i in 1:simd_lane_width:N_SIMD
        @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
        f!(lane+i, v, work, args...)
      end
    end
    # Do the remainder
    for i in N_SIMD+1:N_particle
      @assert last(i) <= N_particle "Out of bounds!"
      f!(i, v, work, args...)
    end
  else
    if N_particle >= multithread_threshold
      Threads.@threads for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    else
      @simd for i in 1:N_particle
        @assert last(i) <= N_particle "Out of bounds!"
        f!(i, v, work, args...)
      end
    end
  end
  return v
end


@generated function genlaunch!(
  chain::C,
  v::A,
  work;
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {C<:NTuple,A}

  # Manually unroll the chain before
  N_fc = length(chain.parameters)
  inner_loop = [
    quote 
      f! = chain[$k][1]
      args = Base.tail(chain[$k])
      f!(i, v, work, args...)
    end for k in 1:N_fc
  ]

  return quote
    N_particle = size(v, 1)

    if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
      lane = VecRange{simd_lane_width}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        for j in 1:simd_lane_width:N_SIMD
          i = lane+j
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          $(inner_loop...)
        end
      else
        for j in 1:simd_lane_width:N_SIMD
          i = lane+j
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          $(inner_loop...)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!" 
        $(inner_loop...)
      end
    else
      if N_particle >= multithread_threshold
        for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          $(inner_loop...)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          $(inner_loop...)
        end
      end
    end
    return v
  end
end

@generated function genlaunch!(
  chain::C,
  v::A,
  work;
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {C<:NTuple,A}

  # Manually unroll the chain before
  N_fc = length(chain.parameters)
  inner_loop = [
    quote 
      f! = chain[$k][1]
      args = Base.tail(chain[$k])
      f!(i, v, work, args...)
    end for k in 1:N_fc
  ]

  return quote
    N_particle = size(v, 1)

    if A <: SIMD.FastContiguousArray && eltype(A) <: SIMD.ScalarTypes && simd_lane_width != 0 # do SIMD
      lane = VecRange{simd_lane_width}(0)
      rmn = rem(N_particle, simd_lane_width)
      N_SIMD = N_particle - rmn
      if N_particle >= multithread_threshold
        for j in 1:simd_lane_width:N_SIMD
          i = lane+j
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          $(inner_loop...)
        end
      else
        for j in 1:simd_lane_width:N_SIMD
          i = lane+j
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          $(inner_loop...)
        end
      end
      # Do the remainder
      for i in N_SIMD+1:N_particle
        @assert last(i) <= N_particle "Out of bounds!" 
        $(inner_loop...)
      end
    else
      if N_particle >= multithread_threshold
        for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          $(inner_loop...)
        end
      else
        @simd for i in 1:N_particle
          @assert last(i) <= N_particle "Out of bounds!"
          $(inner_loop...)
        end
      end
    end
    return v
  end
end