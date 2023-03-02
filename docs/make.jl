push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, Rayons

makedocs(
    modules = [Rayons],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Laurent Plagne",
    sitename = "Rayons.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/LaurentPlagne/Rayons.jl.git",
)
