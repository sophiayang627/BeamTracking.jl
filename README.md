# BeamTracking

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://bmad-sim.github.io/BeamTracking.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://bmad-sim.github.io/BeamTracking.jl/dev/)
[![Build Status](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmad-sim/BeamTracking.jl/actions/workflows/CI.yml?query=branch%3Amain)


To develop this package, you will need [`AcceleratorLattice.jl`](https://github.com/bmad-sim/AcceleratorLattice.jl) and [`AtomicAndPhysicalConstants.jl`](https://github.com/bmad-sim/AtomicAndPhysicalConstants.jl), which is not yet on the official Julia registry. Therefore in the Julia REPL run:

```julia
import Pkg;
Pkg.add(url="https://github.com/bmad-sim/AcceleratorLattice.jl.git");     # Dependency
Pkg.add(url="https://github.com/bmad-sim/AtomicAndPhysicalConstants.jl.git"); # Dependency
Pkg.develop(url="https://github.com/bmad-sim/BeamTracking.jl.git");       # This package! Replace bmad-sim with your username if working on a fork
```

If working on your own fork, replace `bmad-sim` in the above `develop` url with your Github username.

In your `~/.julia/dev/` directory, you will now see the directories `AcceleratorLattice` and `BeamTracking`. Both of these are Github repos where you can do your work and push changes.

If you would like to add a dependency to the project (e.g. [`OrdinaryDiffEq.jl`](https://github.com/SciML/OrdinaryDiffEq.jl)), you should:

1. In the `~/.julia/dev/BeamTracking.jl` directory, start `julia`
2. Enter package-mode using `]`, and then type `activate .`. This activates the `Project.toml` in the current directory you're in
3. In package-mode, type `add OrdinaryDiffEq`. This will add the package as a dependency to the `Project.toml`.
4. In the main module `src/BeamTracking.jl`, add a `using OrdinaryDiffEq` to the top of `BeamTracking.jl`, and in your particularly module defined below add a `using ..OrdinaryDiffEq` which basically says, go one module level up and get `OrdinaryDiffEq`.