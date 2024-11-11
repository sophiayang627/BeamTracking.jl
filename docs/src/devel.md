# Developer's Guide

The design of `BeamTracking.jl` is built with modularity, high performance, differentiability, and polymorphism as key principles. Adding a new type of tracking method or lattice element should be simple.

The entire package is centered around one single `track!` function, which has the following format:

```julia
track!(beam::Beam, element, work=work)
```

Here, `beam` is a `Beam` struct, which is described in detail [below](@ref beam). `element` is some element which to track the beam through, and `work` is an optional tuple of the minimal number of temporaries needed for use inside the tracking function (for some elements, it is an empty tuple). 

After calling `track!`, the `beam` struct is mutated to contain the particle phase space coordinates (and possiblly spin transport quaternions) after propagation through the `element`. With this `track!` function, all one needs to do is define their own custom `element` type, and then when looping through a vector of elements, Julia's multiple dispatch will take care of each particular element.

For example, using the `Linear.Drift` and `Linear.Quadrupole` elements:

```julia
d = Linear.Drift(L=0.2)
qf = Linear.Quadrupole(B1=-18.5, L=0.5)
qd = Linear.Quadrupole(B1=18.5, L=0.5)

fodo = (qf, d, qd, d)

beam = Beam(x=1e-3, px=1e-4, pz=1e-4, beta_gamma_ref=35000.) # Creates a Beam with one particle
for ele in fodo
  track!(beam, ele)
end
```

## [The `Beam` Struct](@id beam)

The `Beam` struct contains the following fields:

- **`species::Species`** -- A `Species` type storing the beam species (e.g. electron)
- **`beta_gamma_ref::Float64`**   -- A reference Lorentz $\beta\gamma$ for normalizing the transverse momenta to
- **`v::T`** -- All particles' 6D phase space coordinates stored in a structure-of-arrays (SoA) memory layout
- **`q::U`** -- If spin tracking, then all particles' spin transport quaternions (using the quaternion type defined in [`ReferenceFrameRotations.jl`](https://github.com/JuliaSpace/ReferenceFrameRotations.jl)) stored in a structure-of-arrays layout. Else, `nothing`

If you are unfamiliar with structure-of-arrays (SoA) and array-of-structures (AoS), you should read [this Wikipedia article](https://en.wikipedia.org/wiki/AoS_and_SoA).

This package extensively uses the [`StructArrays.jl`](https://github.com/JuliaArrays/StructArrays.jl), which features efficient and highly convenient implementation of structs for an SoA memory layout. A basic understanding of this package would be useful too.

In `BeamTracking.jl`, we define the `Coord` type which is a simple [*static* vector](https://github.com/JuliaArrays/StaticArrays.jl) of the 6D phase space coordinates:

```julia
# Static phase space coordinate vector
Base.@kwdef struct Coord{T} <: FieldVector{6, T} 
  x::T  = 0.0
  px::T = 0.0
  y::T  = 0.0
  py::T = 0.0
  z::T  = 0.0
  pz::T = 0.0
end
```

Note that this is *static* and immutable, so you cannot change a field in `Coord`. Also note that a `Vector{Coord{T}}` is an AoS memory layout, and NOT SoA! `StructArrays.jl` makes it easy for us to use this type in a SoA layout, by simple using the `StructArray` type. For example,

```julia
N_particles = 1000
x  = rand(N_particles)
px = rand(N_particles)
y  = rand(N_particles)
py = rand(N_particles)
z  = rand(N_particles)
pz = rand(N_particles)

v = StructArray{Coord{Float64}}((x, px, y, py, z, pz))
v.x  # accesses x array
v.px # accesses px array
v[1] # This goes from SoA to a single Coord struct representing the first Coord! 
```
In the above example, `v` is what is stored in the `Beam` struct. The choice of SoA was made after careful speed benchmarks and the desire to have a mutable interface (`track!` instead of `track.` for a beam).

When there spin tracking, `q`  is a `StructArray{<:Quaternion}`

## Rules for `track!` Implementations

Vectorized operations should be used. The function should be type-stable and, when pre-allocating the necessary `work`, have zero allocations included when tracking a single _non-parametric_ GTPSA map or tracking a beam of particles as immutable numbers (`Float64` or `Dual` numbers for example). For parametric GTPSA maps (e.g. when a quadrupole strength is included as a parameter in the GTPSA), the allocation restriction is loosened in order to maintain readable code. Tests for this are included in the `test_matrix` function in [`runtests.jl`](https://github.com/bmad-sim/BeamTracking.jl/blob/main/test/runtests.jl)

## Compatibility with `GTPSA.jl`

[`GTPSA.jl`](https://github.com/bmad-sim/GTPSA.jl) is a package which is used extensively in the `SciBmad` ecosystem for Taylor map tracking and normal form calculations. 

## Adding a dependency to a submodule

If you would like to add a dependency to the project (e.g. [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl)), you should:

1. In the `~/.julia/dev/BeamTracking.jl` directory, start `julia`
2. Enter package-mode using `]`, and then type `activate .`. This activates the `Project.toml` in the current directory you're in
3. In package-mode, type `add OrdinaryDiffEq`. This will add the package as a dependency to the `Project.toml`.
4. In the main module `src/BeamTracking.jl`, add a `using OrdinaryDiffEq` to the top of `BeamTracking.jl`, and in your particularly module defined below add a `using ..OrdinaryDiffEq` which basically says, go one module level up and get `OrdinaryDiffEq`.
