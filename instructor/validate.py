"""Validate reference and starter packages in disposable directories.

Run with the required Julia packages available:
    python3 instructor/validate.py
Logs are retained in instructor/validation-output. The working-tree starter
and reference solution are never replaced. No files are published.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile
from build_student_package import ROOT, student_files

LOGS = ROOT / "instructor/validation-output"
LOGS.mkdir(exist_ok=True)
REPORTS = {f"{stem}.{ext}" for stem in ("production-versus-budget", "supply-ranges", "range-width-versus-budget")
           for ext in ("csv", "png", "svg")} | {"amino-acid-allocation.csv"}

def stage(reference=False):
    """Copy the student package into a new temporary directory and return its Path.

    If reference is True, replace only that copy's Compute.jl and responses.md
    with the reference files. Keep the temporary directory for inspection and
    leave the working tree unchanged. Copying errors propagate to the caller.
    """
    target = Path(tempfile.mkdtemp(prefix="ps3-validation-"))
    for path in student_files():
        dest = target / path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    if reference:
        for name in ("src/Compute.jl", "responses.md"):
            shutil.copy2(ROOT / "solution" / name, target / name)
    return target

def set_track(target, track):
    """Write track ("standard", "advanced", or an invalid word) to target/TRACK.txt."""
    (target / "TRACK.txt").write_text(track + "\n")

def run(target, label, expected, script="check_submission.jl", args=()):
    """Run a Julia script in target and return its combined stdout/stderr as text.

    target is the staged assignment Path, label names its log, expected is the
    required process exit code, and script is relative to target. Append args
    after the script filename (for example, ("--solution",)). Inherit the
    current environment, write LOGS/<label>.log, and raise AssertionError if the
    exit code differs. Subprocess launch and file errors propagate.
    """
    result = subprocess.run(["julia", "--startup-file=no", "--project=.", script, *args],
                            cwd=target, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, env=os.environ.copy())
    (LOGS / (label + ".log")).write_text(result.stdout)
    assert result.returncode == expected, f"{label}: exit {result.returncode}; see its log"
    print(f"{label}: expected exit {expected}", flush=True)
    return result.stdout

run(ROOT, "solver-interface", 0, script="instructor/test_solver.jl")
# The direct instructor mode must preserve student files and their output files.
protected = [ROOT / "src/Compute.jl", ROOT / "responses.md", ROOT / "TRACK.txt", ROOT / "MANIFEST.txt"]
protected += [p for p in (ROOT / "outputs").rglob("*") if p.is_file()]
before = {p: p.read_bytes() if p.exists() else None for p in protected}
# The reference solution checks the Advanced track even though TRACK.txt says standard.
output = run(ROOT, "solution-flag-checker", 0, args=("--solution",))
assert output.count("24     24") == 3, output
assert "REFERENCE SOLUTION PASSED" in output and "READY TO PACKAGE" not in output
manifest = (ROOT / "solution/MANIFEST.txt").read_text()
assert "solution/src/Compute.jl" in manifest and "solution/responses.md" in manifest
output = run(ROOT, "solution-flag-part-1", 0, script="testme_part_1.jl", args=("--solution",))
assert "24     24" in output
run(ROOT, "solution-flag-driver", 0, script="runproduction.jl", args=("--solution",))
assert (ROOT / "solution/outputs/supply-ranges.csv").is_file()
assert all((p.read_bytes() if p.exists() else None) == content for p, content in before.items())

reference = stage(reference=True)
# A byte-order mark, capital letter, and Windows line ending must still select Advanced.
(reference / "TRACK.txt").write_bytes(b"\xef\xbb\xbfAdvanced\r\n")
output = run(reference, "reference-checker", 0)
assert output.count("24     24") == 3, output
assert "READY TO PACKAGE" in output
preview = ROOT / "instructor/reference"
preview.mkdir(exist_ok=True)
for path in (reference / "outputs").iterdir():
    shutil.copy2(path, preview / path.name)
assert "Track checked: advanced" in (reference / "MANIFEST.txt").read_text()
# The checker clears exactly these reports, so the list must match what it writes.
assert {p.name for p in (reference / "outputs").iterdir()} == REPORTS
set_track(reference, "standard")

starter = stage()
for script in ("check_submission.jl", "testme_part_1.jl", "runproduction.jl"):
    output = run(starter, "missing-solution-" + script, 1, script=script, args=("--solution",))
    assert "--solution requires" in output
assert not (starter / "outputs").exists() and not (starter / "MANIFEST.txt").exists()
output = run(starter, "starter-checker", 1)
assert "READY TO PACKAGE" not in output
assert all(f"Question {i} still has a placeholder" in output for i in (1, 2))
assert "Question 3 still has a placeholder" not in output
assert "Running testme_part_3.jl" not in output
assert (starter / "MANIFEST.txt").is_file()
assert not (starter / "outputs").exists()
# An invalid track must stop the checker and driver with a clear message.
set_track(starter, "expert")
for script in ("check_submission.jl", "runproduction.jl"):
    output = run(starter, "invalid-track-" + script, 1, script=script)
    assert "TRACK.txt must contain standard or advanced" in output

# A partially completed submission should keep its Part 1 credit.
partial = stage()
source = (ROOT / "solution/src/Compute.jl").read_text()
start = source.index("function production_curve")
source = source[:start] + 'function production_curve(model, budgets::AbstractVector{<:Real})\n    error("unfinished curve");\nend\n'
(partial / "src/Compute.jl").write_text(source)
# Earlier reports must not survive a failing run; other files in outputs/ stay.
(partial / "outputs").mkdir()
for name in ("production-versus-budget.png", "supply-ranges.csv", "notes.txt"):
    (partial / "outputs" / name).write_text("stale")
output = run(partial, "partial-checker", 1)
assert [p.name for p in (partial / "outputs").iterdir()] == ["notes.txt"]
assert "testme_part_1.jl: all tests passed" in output
assert "testme_part_2.jl: some tests failed" in output
assert "Running testme_part_3.jl" not in output

# Completed required work passes with the Advanced function and answer untouched.
# Selecting the Advanced track with unfinished work must preserve the required-work status.
parts12 = stage()
source = (ROOT / "solution/src/Compute.jl").read_text()
start = source.index('\n"""\n    supply_ranges(')
stub = (ROOT / "src/Compute.jl").read_text()
source = source[:start] + stub[stub.index('\n"""\n    supply_ranges('):]
(parts12 / "src/Compute.jl").write_text(source)
answers = (ROOT / "solution/responses.md").read_text()
advanced_prompt = (ROOT / "responses.md").read_text().split('3. **', 1)[1]
(parts12 / "responses.md").write_text(answers.split('3. **', 1)[0] + '3. **' + advanced_prompt)
output = run(parts12, "parts-1-2-checker", 0)
assert "testme_part_1.jl: all tests passed" in output
assert "testme_part_2.jl: all tests passed" in output
assert "Running testme_part_3.jl" not in output
assert "Question 3 still has a placeholder" not in output
assert "READY TO PACKAGE" in output
assert "Advanced: not selected in TRACK.txt (optional)" in (parts12 / "MANIFEST.txt").read_text()
assert (parts12 / "outputs/production-versus-budget.png").is_file()
assert not (parts12 / "outputs/supply-ranges.csv").exists()
output = run(parts12, "parts-1-2-driver", 0, script="runproduction.jl")
assert "Standard track selected in TRACK.txt" in output
assert (parts12 / "outputs/production-versus-budget.png").is_file()
assert not (parts12 / "outputs/supply-ranges.csv").exists()
set_track(parts12, "advanced")
output = run(parts12, "unfinished-advanced-checker", 1)
assert "testme_part_3.jl: some tests failed" in output
assert "Required work: checks passed" in output
assert "Advanced discussion: Question 3 still has a placeholder" in output
assert "one Magic Point pending" not in output
output = run(parts12, "unfinished-advanced-driver", 1, script="runproduction.jl")
assert "Complete supply_ranges" in output

# Correct Advanced code alone cannot satisfy the Advanced writing requirement.
shutil.copy2(ROOT / "solution/src/Compute.jl", parts12 / "src/Compute.jl")
output = run(parts12, "missing-advanced-answer", 1)
assert "testme_part_3.jl: all tests passed" in output
assert "Required work: checks passed" in output
assert "Advanced discussion: Question 3 still has a placeholder" in output
assert "one Magic Point pending" not in output

# An Advanced report failure must not invalidate the required reports or status.
shutil.copy2(ROOT / "solution/responses.md", parts12 / "responses.md")
reports = parts12 / "src/Reporting.jl"
with reports.open('a') as handle:
    handle.write('\n@eval ProductionReports function write_range_outputs(assignment, output_directory::AbstractString)\n'
                 '    error("Deliberate Advanced reporting failure")\nend\n')
output = run(parts12, "advanced-reporting-error", 1)
assert "Required work: checks passed" in output
assert "Advanced outputs: OUTPUT GENERATION NEEDS ATTENTION" in output
assert "one Magic Point pending" not in output

# Load errors must remain distinct from test failures and still yield a manifest.
syntax = stage()
(syntax / "src/Compute.jl").write_text("function unfinished(\n")
output = run(syntax, "syntax-error-checker", 1)
assert "tests could not run" in output
assert (syntax / "MANIFEST.txt").is_file()

# Correct code with missing writing must not be reported ready.
(reference / "responses.md").unlink()
output = run(reference, "missing-responses-checker", 1)
assert "responses.md is missing" in output and "READY TO PACKAGE" not in output

# A rendering failure must not erase passing test results or claim readiness.
shutil.copy2(ROOT / "solution/responses.md", reference / "responses.md")
(reference / "src/Reporting.jl").write_text('error("Deliberate reporting failure")\n')
output = run(reference, "reporting-error-checker", 1)
assert "OUTPUT GENERATION NEEDS ATTENTION" in output
assert "testme_part_1.jl: all tests passed" in output
assert "READY TO PACKAGE" not in output
print("Validated both tracks, invalid tracks, unfinished optional work, syntax errors, missing responses, and reporting failures.")
