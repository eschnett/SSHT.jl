# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using SSHT

makedocs(; sitename="SSHT", format=Documenter.HTML(), modules=[SSHT])

deploydocs(; repo="github.com/eschnett/SSHT.jl.git", devbranch="main", push_preview=true)
