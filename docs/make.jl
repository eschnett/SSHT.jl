# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using FastSphericalHarmonics

makedocs(; sitename="ssht", format=Documenter.HTML(), modules=[ssht])

deploydocs(; repo="github.com/eschnett/ssht.jl.git", devbranch="main", push_preview=true)
