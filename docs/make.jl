using BeamTracking
using Documenter

makedocs(;
    authors="mattsignorelli <mgs255@cornell.edu> and contributors",
    sitename="BeamTracking.jl",
  format=Documenter.HTMLWriter.HTML(size_threshold = nothing),
  pages = 
  [
    "Home" => "index.md",
    "Submodules" => ["linear.md",
                     "matrixkick.md",
                     "misc.md"],
    "Developer's Guide" => "devel.md"
  ]
)

deploydocs(;
    repo="github.com/bmad-sim/BeamTracking.jl.git",
)
