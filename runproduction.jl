# ----------------------------------------------------------------------------
# Generate the production and supply-range reports
# ----------------------------------------------------------------------------
# Run from the assignment root:
#   julia --startup-file=no --project=. runproduction.jl
#   julia --startup-file=no --project=. runproduction.jl --solution
#
# The --solution flag loads the reference code, writes every report, and writes
# solution/outputs/.
#
# After the Part 2 functions are complete, write the production plot and data
# tables to outputs/. When TRACK.txt says advanced, also write the optional Part 3
# supply-range plots and tables. Existing files with the same names are overwritten.


# 1. Load the implementation and reporting utilities -
ENV["GKSwstype"] = "100"; # render plots to files without opening a GUI window

include(joinpath(@__DIR__, "Include.jl"));
include(joinpath(@__DIR__, "src", "Reporting.jl"));


# 2. Select the track and output directory -
track_path = joinpath(@__DIR__, "TRACK.txt");
isfile(track_path) || error("TRACK.txt is missing; restore it with the word standard or advanced.");
track = lowercase(strip(replace(read(track_path, String), '\ufeff' => "")));  # drop a Windows byte-order mark
track in ("standard", "advanced") || error("TRACK.txt must contain standard or advanced, not $(repr(track)).");
advanced_mode = "--solution" in ARGS || track == "advanced";

output_directory = "--solution" in ARGS ?
    joinpath(@__DIR__, "solution", "outputs") : joinpath(@__DIR__, "outputs");

"--solution" in ARGS && println("Running the instructor reference solution.");


# 3. Write the Part 2 production plot and allocation tables -
ProductionReports.write_outputs(CellFreeProduction, output_directory);


# 4. Write the optional Advanced supply-range plots and tables -
if advanced_mode
    ProductionReports.write_range_outputs(CellFreeProduction, output_directory);
else
    println("Standard track selected in TRACK.txt; skipping the Part 3 supply-range reports.");
end
