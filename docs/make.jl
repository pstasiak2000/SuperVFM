# push!(LOAD_PATH,"../src/")
using Documenter
using SuperVFM

#Generates the subpages of the vortex initial conditions
SUBSECTION_PAGES = "vortex_configs/" .* readdir("./docs/src/vortex_configs/")

makedocs(sitename="SuperVFM.jl", remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Vortex Configurations" => SUBSECTION_PAGES
            ]    
)

deploydocs(
    repo = "github.com/pstasiak2000/SuperVFM.jl.git",
)

