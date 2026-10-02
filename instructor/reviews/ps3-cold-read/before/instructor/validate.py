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

def stage(reference=False):
    target = Path(tempfile.mkdtemp(prefix="ps3-validation-"))
    for path in student_files():
        dest = target / path.relative_to(ROOT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    if reference:
        for name in ("src/Compute.jl", "responses.md"):
            shutil.copy2(ROOT / "solution" / name, target / name)
    return target

def run(target, label, expected, script="check_submission.jl"):
    result = subprocess.run(["julia", "--startup-file=no", "--project=.", script],
                            cwd=target, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, env=os.environ.copy())
    (LOGS / (label + ".log")).write_text(result.stdout)
    assert result.returncode == expected, f"{label}: exit {result.returncode}; see its log"
    print(f"{label}: expected exit {expected}", flush=True)
    return result.stdout

reference = stage(reference=True)
output = run(reference, "reference-checker", 0)
assert output.count("24     24") == 2, output
assert "READY TO PACKAGE" in output
preview = ROOT / "instructor/reference"
preview.mkdir(exist_ok=True)
for path in (reference / "outputs").iterdir():
    shutil.copy2(path, preview / path.name)
assert (reference / "MANIFEST.txt").is_file()

starter = stage()
output = run(starter, "starter-checker", 1)
assert "READY TO PACKAGE" not in output
assert all(f"Question {i} still has a placeholder" in output for i in range(1, 4))
assert (starter / "MANIFEST.txt").is_file()
assert not (starter / "outputs").exists()

# A partially completed submission should keep its Part 1 credit.
partial = stage()
source = (ROOT / "solution/src/Compute.jl").read_text()
start = source.index("function production_curve")
source = source[:start] + 'function production_curve(model, budgets::AbstractVector{<:Real})\n    error("unfinished curve");\nend\n'
(partial / "src/Compute.jl").write_text(source)
output = run(partial, "partial-checker", 1)
assert "testme_part_1.jl: all tests passed" in output
assert "testme_part_2.jl: some tests failed" in output

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
print("Validated reference, starter, partial implementation, syntax errors, missing responses, and reporting failures.")
