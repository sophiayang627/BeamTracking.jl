
const REGISTER_SIZE = VectorizationBase.register_size()
const XI  = 1
const PXI = 2
const YI  = 3
const PYI = 4
const ZI  = 5
const PZI = 6

#=
@generated function foo(chain::T) where {T<:Tuple}
  N = length(T.parameters)
  exprs = []
  for i in 1:N
      push!(exprs, quote
          f = chain[$i][1]
          args = Base.tail(chain[$i])  # same as chain[$i][2:end] but type-stable
          println(f, args...)
      end)
  end
  return Expr(:block, exprs...)
end

@generated function foo(chain::T) where {T<:Tuple}
  N = length(T.parameters)

  # Any other setup outside the loop

  # Build loop body
  loop_body = [
      quote
          f = chain[$i][1]
          args = Base.tail(chain[$i])
          println("iteration $i")
          f(args...)
      end
      for i in 1:N
  ]

  return quote
      for i in 1:10
        $(loop_body...)
      end
      println("Done.")
  end
end
=#

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
function launch!(
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





#=
function launch!(
  chain::C,
  v::A,
  work;
  simd_lane_width=0, # autovectorize by default #floor(Int, REGISTER_SIZE/sizeof(eltype(A))),
  multithread_threshold=Threads.nthreads() > 1 ? 1750*Threads.nthreads() : typemax(Int),
) where {C<:Tuple,A}

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
        Threads.@threads for j in 1:simd_lane_width:N_SIMD
          i = lane+j
          @assert last(i) <= N_particle "Out of bounds!"  # Use last because VecRange SIMD
          $(inner_loop...)
        end
      else
        Threads.@threads for j in 1:simd_lane_width:N_SIMD
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
        Threads.@threads for i in 1:N_particle
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
=#