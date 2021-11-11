using Documenter
using OneDmin

makedocs(
    sitename = "OneDmin.jl",
    format = Documenter.HTML(assets = ["assets/style.css"], prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [OneDmin],
    pages = ["Home" => "index.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/vepiteski/OneDmin.jl")
#https://juliadocs.github.io/Documenter.jl/stable/man/hosting/ ?
