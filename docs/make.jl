push!(LOAD_PATH, "../src/")
using Documenter, EigenstratFormat

makedocs(
    sitename = "EigenstratFormat.jl",
    modules = [EigenstratFormat],
#    remotes = nothing
)
deploydocs(;
    repo="github.com/yogischogi/EigenstratFormat.jl",
    devbranch = "main"
)


