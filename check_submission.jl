# Run all three test suites, check the responses, generate outputs, and write a manifest.
# Usage: julia --startup-file=no --project=. check_submission.jl
# Exit with code 0 when the tests, output generation, and response checks pass;
# otherwise exit with code 1. A score still requires instructor review.
using Test;  # distinguish failed assertions from source-loading errors
using SHA;   # fingerprint submitted source files and responses
using Dates; # timestamp the submission manifest

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
    discussion_issues(path) -> Vector{String}

Read the response file at `path` and return messages for missing or duplicate
questions, empty answers, and remaining TODO/TBD placeholders. Expect the five
numbered, bold question prompts followed by answer paragraphs in responses.md.
Ignore HTML comments when checking whether an answer is present.

Return an empty vector if no issues are detected; this does not grade correctness.
A missing or unreadable file returns a warning message. The file is unchanged.
"""
function discussion_issues(path)
    isfile(path) || return ["responses.md is missing"];
    content = try
        replace(read(path, String), "\r\n" => "\n");
    catch
        return ["responses.md could not be read"];
    end
    answers = Dict{Int,String}();
    issues = String[];
    # Separate each numbered prompt from the answer paragraphs that follow it -
    for section in eachmatch(r"(?ms)^[ \t]{0,2}([1-5])\.[ \t]+\*\*.*?(?=^[ \t]{0,2}[1-5]\.[ \t]+\*\*|\z)", content)
        number = parse(Int, section.captures[1]);
        haskey(answers, number) && push!(issues, "Question $number appears more than once");
        paragraphs = split(section.match, r"\n[ \t]*\n"; limit=2);
        answer = length(paragraphs) == 2 ? paragraphs[2] : "";
        answers[number] = strip(replace(answer, r"(?s)<!--.*?-->" => ""));
    end
    for number in 1:5
        answer = get(answers, number, "");
        if !occursin(r"[\p{L}\p{N}]", answer)
            push!(issues, "Question $number has no answer");
        elseif occursin(r"(?i)\b(TODO|TBD)\b|Write your response|Your answer here", answer)
            push!(issues, "Question $number still has a placeholder");
        end
    end
    return issues;
end

println("PS3 submission check: 24 tests per part, 72 total.");
statuses = Dict{String,String}();
suites = Dict{String,Module}();
# Isolate each suite so loading the same student module twice does not replace it -
for filename in ("testme_part_1.jl", "testme_part_2.jl", "testme_part_3.jl")
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

issues = discussion_issues(joinpath(@__DIR__, "responses.md"));
all_passed = all(==("all tests passed"), values(statuses));
report_status = "not generated because tests need attention";
if all_passed
    try
        ENV["GKSwstype"] = "100";
        include(joinpath(@__DIR__, "src", "Reporting.jl"));
        assignment = getfield(suites["testme_part_3.jl"], :CellFreeProduction);
        ProductionReports.write_outputs(assignment, joinpath(@__DIR__, "outputs"));
        ProductionReports.write_range_outputs(assignment, joinpath(@__DIR__, "outputs"));
        global report_status = "plot and data tables generated";
    catch caught
        global report_status = "OUTPUT GENERATION NEEDS ATTENTION";
        showerror(stderr, caught, catch_backtrace());
        println(stderr);
    end
end

# Record the test status and fingerprint every submitted source file -
files = [joinpath(@__DIR__, "Include.jl"), joinpath(@__DIR__, "responses.md")];
for (directory, _, filenames) in walkdir(joinpath(@__DIR__, "src"))
    append!(files, [joinpath(directory, name) for name in filenames]);
end
open(joinpath(@__DIR__, "MANIFEST.txt"), "w") do io
    println(io, "PS3 CHEME 4800/5800 Fall 2026 submission manifest");
    println(io, "Generated: ", Dates.now());
    for name in sort(collect(keys(statuses)))
        println(io, name, ": ", statuses[name]);
    end
    println(io, "Outputs: ", report_status);
    println(io, "Discussion: ", isempty(issues) ? "no unanswered prompts detected; instructor review required" : join(issues, "; "));
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
ready = all_passed && report_status == "plot and data tables generated" && isempty(issues);
println("Status: ", ready ? "READY TO PACKAGE; instructor review pending" : "WORK NEEDS ATTENTION");
println("Wrote MANIFEST.txt. Zip the assignment folder, including responses and outputs.");
println("Name it CHEME-4800-5800-PS3-<your netid>.zip and upload it to Canvas yourself.");
println("Submit your current work by the initial deadline even if some tests still fail.");
exit(ready ? 0 : 1);
