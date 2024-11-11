# Developer's Guide

The design of `BeamTracking.jl` is built with modularity, high performance, differentiability, and polymorphism as key principles. Adding a new type of tracking method or lattice element should be simple.

The entire package is centered around one single `track!` function, which has the following format:

```julia
track!(beam::Beam, element, work=work)
```

Here, `beam` is a `Beam` struct, which is described in detail [below](@ref Beam). `element` is some element which to track the beam through, and `work` is an optional tuple of the minimal number of temporaries needed for use inside the tracking function (for some elements, it is an empty tuple). 

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

## [The `Beam` Struct](@id Beam)

The `Beam` struct contains



## Adding a dependency to a submodule