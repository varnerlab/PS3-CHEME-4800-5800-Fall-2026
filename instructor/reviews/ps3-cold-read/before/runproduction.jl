# Run from the assignment root: julia --startup-file=no --project=. runproduction.jl
ENV["GKSwstype"] = "100"; # render plots to files without opening a GUI window
include(joinpath(@__DIR__, "Include.jl"));
include(joinpath(@__DIR__, "src", "Reporting.jl"));
ProductionReports.write_outputs(CellFreeProduction, joinpath(@__DIR__, "outputs"));
