# Check required work, generate outputs, and write a manifest.
# Usage: julia --startup-file=no --project=. check_submission.jl
# TRACK.txt selects standard or advanced; advanced also checks Part 3 and question 3.
# Add --solution to check local reference code and responses on the Advanced track;
# write its outputs and manifest under solution/. Student source files and responses stay unchanged.
# Exit with code 0 when the tests, output generation, and response checks pass;
# otherwise exit with code 1. A score still requires instructor review.
using Test;  # distinguish failed assertions from source-loading errors
using SHA;   # fingerprint submitted source files and responses
using Dates; # timestamp the submission manifest

# Choose the files being checked before loading any suite or writing outputs -
const SOLUTION_MODE = "--solution" in ARGS;
# Read the selected track; the reference solution always checks the Advanced track -
const TRACK_PATH = joinpath(@__DIR__, "TRACK.txt");
isfile(TRACK_PATH) || error("TRACK.txt is missing; restore it with the word standard or advanced.");
const TRACK = lowercase(strip(replace(read(TRACK_PATH, String), '\ufeff' => "")));  # drop a Windows byte-order mark
TRACK in ("standard", "advanced") || error("TRACK.txt must contain standard or advanced, not $(repr(TRACK)).");
const ADVANCED_MODE = SOLUTION_MODE || TRACK == "advanced";
const CHECK_ROOT = SOLUTION_MODE ? joinpath(@__DIR__, "solution") : @__DIR__;
const RESPONSE_PATH = joinpath(CHECK_ROOT, "responses.md");
const COMPUTE_PATH = joinpath(CHECK_ROOT, "src", "Compute.jl");
if SOLUTION_MODE
    for path in (COMPUTE_PATH, RESPONSE_PATH)
        isfile(path) || error("--solution requires the local instructor file: $path");
    end
    println("Checking the instructor reference solution; outputs and manifest go under solution/.");
end

"""
    rootcause(caught)

Inspect a caught exception and return the first exception that is not a
`LoadError`, unwrapping nested file-include errors. Return `caught` itself if
it has no wrapper. Leave the exception unchanged. The checker uses this result
to distinguish failed test assertions from errors that prevented a suite running.
"""
function rootcause(caught)
    while caught isa LoadError
        caught = caught.error;
    end
    return caught;
end

"""
    discussion_issues(path; questions=1:2) -> Vector{String}

Read the response file at `path` and return messages for missing or duplicate
questions, empty answers, and remaining TODO/TBD placeholders for `questions`.
Expect numbered, bold question prompts followed by answer paragraphs in responses.md.
Ignore HTML comments when checking whether an answer is present.

Return an empty vector if no issues are detected; this does not grade correctness.
A missing or unreadable file returns a warning message. The file is unchanged.
"""
function discussion_issues(path; questions=1:2)
    isfile(path) || return ["responses.md is missing"];
    content = try
        replace(read(path, String), "\r\n" => "\n");
    catch
        return ["responses.md could not be read"];
    end
    answers = Dict{Int,String}();
    issues = String[];
    # Separate each numbered prompt from the answer paragraphs that follow it -
    for section in eachmatch(r"(?ms)^[ \t]{0,2}([1-3])\.[ \t]+\*\*.*?(?=^[ \t]{0,2}[1-3]\.[ \t]+\*\*|\z)", content)
        number = parse(Int, section.captures[1]);
        number in questions || continue;
        haskey(answers, number) && push!(issues, "Question $number appears more than once");
        paragraphs = split(section.match, r"\n[ \t]*\n"; limit=2);
        answer = length(paragraphs) == 2 ? paragraphs[2] : "";
        answers[number] = strip(replace(answer, r"(?s)<!--.*?-->" => ""));
    end
    for number in questions
        answer = get(answers, number, "");
        if !occursin(r"[\p{L}\p{N}]", answer)
            push!(issues, "Question $number has no answer");
        elseif occursin(r"(?i)\b(TODO|TBD)\b|Write your response|Your answer here", answer)
            push!(issues, "Question $number still has a placeholder");
        end
    end
    return issues;
end

println("PS3 submission check, ", ADVANCED_MODE ? "Advanced" : "Standard", " track: Parts 1 and 2, 48 required tests.");
println(ADVANCED_MODE ? "Advanced: also checking Part 3 (24 tests) and question 3." :
    "Advanced: not selected in TRACK.txt; Part 3 and question 3 are optional.");
statuses = Dict{String,String}();
suites = Dict{String,Module}();
required_suites = ("testme_part_1.jl", "testme_part_2.jl");
selected_suites = ADVANCED_MODE ? (required_suites..., "testme_part_3.jl") : required_suites;
# Isolate each suite so loading the same student module twice does not replace it -
for filename in selected_suites
    println("\nRunning ", filename);
    try
        suite = Module(gensym(:PS3Part));
        Base.include(suite, joinpath(@__DIR__, filename));
        suites[filename] = suite;
        statuses[filename] = "all tests passed";
    catch caught
        if rootcause(caught) isa Test.TestSetException
            statuses[filename] = "some tests failed";
        else
            statuses[filename] = "tests could not run";
            showerror(stderr, caught, catch_backtrace());
            println(stderr);
        end
    end
end

issues = discussion_issues(RESPONSE_PATH);
advanced_issues = ADVANCED_MODE ? discussion_issues(RESPONSE_PATH; questions=(3,)) : String[];
required_passed = all(name -> statuses[name] == "all tests passed", required_suites);
advanced_passed = ADVANCED_MODE && statuses["testme_part_3.jl"] == "all tests passed";
# Remove earlier reports so outputs/ holds only results from the code just checked -
for name in ("production-versus-budget.csv", "production-versus-budget.png",
        "production-versus-budget.svg", "amino-acid-allocation.csv", "supply-ranges.csv",
        "supply-ranges.png", "supply-ranges.svg", "range-width-versus-budget.csv",
        "range-width-versus-budget.png", "range-width-versus-budget.svg")
    rm(joinpath(CHECK_ROOT, "outputs", name); force=true);
end
report_status = "not generated because tests need attention";
advanced_report_status = ADVANCED_MODE ? "not generated because tests need attention" : "not requested";
if required_passed
    try
        ENV["GKSwstype"] = "100";
        include(joinpath(@__DIR__, "src", "Reporting.jl"));
        assignment = getfield(suites["testme_part_2.jl"], :CellFreeProduction);
        ProductionReports.write_outputs(assignment, joinpath(CHECK_ROOT, "outputs"));
        global report_status = "plot and data tables generated";
    catch caught
        global report_status = "OUTPUT GENERATION NEEDS ATTENTION";
        showerror(stderr, caught, catch_backtrace());
        println(stderr);
    end
end
if advanced_passed && report_status == "plot and data tables generated"
    try
        assignment = getfield(suites["testme_part_3.jl"], :CellFreeProduction);
        ProductionReports.write_range_outputs(assignment, joinpath(CHECK_ROOT, "outputs"));
        global advanced_report_status = "plot and data tables generated";
    catch caught
        global advanced_report_status = "OUTPUT GENERATION NEEDS ATTENTION";
        showerror(stderr, caught, catch_backtrace());
        println(stderr);
    end
end

required_ready = required_passed && report_status == "plot and data tables generated" && isempty(issues);
advanced_ready = required_ready && advanced_passed &&
    advanced_report_status == "plot and data tables generated" && isempty(advanced_issues);
advanced_status = !ADVANCED_MODE ? "not selected in TRACK.txt (optional)" :
    advanced_ready ? "checks passed; one Magic Point pending teaching-team review" : "WORK NEEDS ATTENTION";

# Record the test status and fingerprint every submitted source file -
files = [joinpath(@__DIR__, "Include.jl"), TRACK_PATH, RESPONSE_PATH];
for (directory, _, filenames) in walkdir(joinpath(@__DIR__, "src"))
    # Fingerprint the implementation actually loaded, including in solution mode -
    append!(files, [joinpath(directory, name) == joinpath(@__DIR__, "src", "Compute.jl") ?
        COMPUTE_PATH : joinpath(directory, name) for name in filenames]);
end
open(joinpath(CHECK_ROOT, "MANIFEST.txt"), "w") do io
    println(io, SOLUTION_MODE ? "PS3 instructor reference solution manifest" :
        "PS3 CHEME 4800/5800 Fall 2026 submission manifest");
    println(io, "Generated: ", Dates.now());
    println(io, "Track checked: ", ADVANCED_MODE ? "advanced" : "standard");
    for name in sort(collect(keys(statuses)))
        println(io, name, ": ", statuses[name]);
    end
    println(io, "Outputs: ", report_status);
    println(io, "Discussion: ", isempty(issues) ? "no unanswered prompts detected; instructor review required" : join(issues, "; "));
    println(io, "Required work: ", required_ready ? "checks passed; instructor review required" : "WORK NEEDS ATTENTION");
    println(io, "Advanced: ", advanced_status);
    if ADVANCED_MODE
        println(io, "Advanced outputs: ", advanced_report_status);
        println(io, "Advanced discussion: ", isempty(advanced_issues) ? "no unanswered prompts detected; instructor review required" : join(advanced_issues, "; "));
    end
    for path in sort(files)
        relative = replace(relpath(path, @__DIR__), '\\' => '/');
        try
            println(io, bytes2hex(open(sha256, path)), "  ", relative);
        catch
            println(io, "MISSING OR UNREADABLE  ", relative);
        end
    end
end

println("\nSubmission summary");
for name in sort(collect(keys(statuses)))
    println(name, ": ", statuses[name]);
end
println("Outputs: ", report_status);
for issue in issues
    println("Discussion: ", issue);
end
println("Required work: ", required_ready ? "checks passed; instructor review pending" : "WORK NEEDS ATTENTION");
println("Advanced: ", advanced_status);
if ADVANCED_MODE
    println("Advanced outputs: ", advanced_report_status);
    for issue in advanced_issues
        println("Advanced discussion: ", issue);
    end
end
ready = required_ready && (!ADVANCED_MODE || advanced_ready);
if SOLUTION_MODE
    println("Status: ", ready ? "REFERENCE SOLUTION PASSED" : "REFERENCE SOLUTION NEEDS ATTENTION");
    println("Wrote solution/MANIFEST.txt; reference outputs are in solution/outputs/.");
else
    println("Status: ", ready ? "READY TO PACKAGE; instructor review pending" : "WORK NEEDS ATTENTION");
    println("Wrote MANIFEST.txt. Zip the assignment folder, including responses and outputs.");
    println("Name it CHEME-4800-5800-PS3-<your netid>.zip and upload it to Canvas yourself.");
    println("Submit your current work by the initial deadline even if some tests still fail.");
end
exit(ready ? 0 : 1);
