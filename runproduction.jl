# Run from the assignment root: julia --startup-file=no --project=. runproduction.jl
# After the Part 2 functions are complete, write the production plot and data
# tables to outputs/. After supply_ranges is complete, also write the Part 3
# supply-range plots and tables. Existing files with the same names are overwritten.
ENV["GKSwstype"] = "100"; # render plots to files without opening a GUI window
include(joinpath(@__DIR__, "Include.jl"));
include(joinpath(@__DIR__, "src", "Reporting.jl"));
ProductionReports.write_outputs(CellFreeProduction, joinpath(@__DIR__, "outputs"));
# Skip only the Part 3 outputs while the supply_ranges starter placeholder remains -
try
    ProductionReports.write_range_outputs(CellFreeProduction, joinpath(@__DIR__, "outputs"));
catch caught
    caught isa ErrorException && occursin("Complete supply_ranges", caught.msg) || rethrow();
    println("Skipped the Part 3 supply-range outputs; complete supply_ranges first.");
end
