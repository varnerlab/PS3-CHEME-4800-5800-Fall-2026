"""Build and test the student solution ZIP using Python's standard library.

Run from the repository root: python3 instructor/build_solution_release.py
Requires the course Julia version (1.12.7) on PATH. The extracted ZIP installs
its recorded Julia environment, so the first run can take several minutes.
"""

import argparse
import hashlib
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_NAME = "PS3-CHEME-4800-5800-Fall-2026-solution"
SUPPORT_FILES = (
    "Include.jl",
    "Project.toml",
    "Manifest.toml",
    "check_submission.jl",
    "runproduction.jl",
    "test_support.jl",
    "testme_part_1.jl",
    "testme_part_2.jl",
    "testme_part_3.jl",
    "RUBRIC.md",
    "LICENSE",
    "THIRD_PARTY_NOTICES.txt",
)
REPORTS = {f"{stem}.{ext}" for stem in ("production-versus-budget", "supply-ranges", "range-width-versus-budget")
           for ext in ("csv", "png", "svg")} | {"amino-acid-allocation.csv"}


def check_markdown(package):
    """Stop the build if a Markdown file links to a missing file or uses math GitHub cannot render."""
    for markdown in package.rglob("*.md"):
        text = markdown.read_text(encoding="utf-8")
        for macro in ("\\operatorname", "\\DeclareMathOperator", "\\newcommand"):
            if macro in text:
                raise ValueError(f"{markdown.name} uses {macro}, which GitHub cannot render")
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if target.startswith(("https:", "http:", "mailto:", "#")):
                continue
            if not (markdown.parent / target.split("#")[0]).exists():
                raise FileNotFoundError(f"Broken local link in {markdown.name}: {target}")


def run_julia(directory, *arguments):
    """Run Julia with the package's environment; fail the build on a nonzero exit."""
    return subprocess.run(
        ["julia", "--startup-file=no", "--project=.", *arguments],
        cwd=directory,
        capture_output=True,
        text=True,
        check=True,
        timeout=3600,
    )


def build_release(output_directory):
    """Return the ZIP path after the public tests and checker pass on its contents."""
    with tempfile.TemporaryDirectory(prefix="ps3-solution-release-") as temporary:
        staging = Path(temporary)
        package = staging / PACKAGE_NAME
        package.mkdir()

        # Package the assignment materials and a README for the completed solution.
        for filename in SUPPORT_FILES:
            shutil.copy2(ROOT / filename, package / filename)
        shutil.copy2(ROOT / "README.md", package / "ASSIGNMENT.md")
        shutil.copy2(ROOT / "instructor" / "solution-README.md", package / "README.md")
        shutil.copytree(ROOT / "data", package / "data", ignore=shutil.ignore_patterns(".DS_Store"))

        # In this package the assignment text is ASSIGNMENT.md, so the rubric links there.
        rubric = (package / "RUBRIC.md").read_text(encoding="utf-8")
        if "[README.md](README.md)" not in rubric:
            raise RuntimeError("RUBRIC.md no longer links to README.md; update the builder")
        rubric = rubric.replace("[README.md](README.md)", "[ASSIGNMENT.md](ASSIGNMENT.md)")
        (package / "RUBRIC.md").write_text(rubric, encoding="utf-8")

        # Put the reference code where Include.jl loads it and keep the starter for comparison.
        # Only Compute.jl differs from the starter; the other source files are supplied code.
        solution = ROOT / "solution" / "src" / "Compute.jl"
        if not solution.is_file():
            raise FileNotFoundError(f"Missing reference solution: {solution}")
        if 'error("Complete ' in solution.read_text(encoding="utf-8"):
            raise RuntimeError("The reference solution still contains a starter placeholder")
        shutil.copytree(ROOT / "src", package / "src", ignore=shutil.ignore_patterns(".DS_Store"))
        (package / "starter").mkdir()
        shutil.copy2(ROOT / "src" / "Compute.jl", package / "starter" / "Compute.jl")
        shutil.copy2(solution, package / "src" / "Compute.jl")

        # The solution covers both tracks, so the checker runs all 72 tests.
        (package / "TRACK.txt").write_text("advanced\n", encoding="utf-8")

        # The model answers move from solution/ to the package root, so their
        # links up one directory become links within the package.
        responses = (ROOT / "solution" / "responses.md").read_text(encoding="utf-8")
        (package / "responses.md").write_text(responses.replace("](../", "]("), encoding="utf-8")
        check_markdown(package)

        # Build one top-level folder so extracting the ZIP keeps its files together.
        candidate = staging / f"{PACKAGE_NAME}.zip"
        with zipfile.ZipFile(candidate, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(package.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(staging))

        # Test the actual archive, including its paths, after installing its environment.
        extracted = staging / "extracted"
        with zipfile.ZipFile(candidate) as archive:
            archive.extractall(extracted)
        student = extracted / PACKAGE_NAME
        run_julia(student, "-e", "using Pkg; Pkg.instantiate()")
        for test_file in ("testme_part_1.jl", "testme_part_2.jl", "testme_part_3.jl"):
            # The student checker catches failures; direct test runs must fail the build.
            result = run_julia(student, test_file)
            print(result.stdout)

        # The checker must accept the model answers and write every report.
        result = run_julia(student, "check_submission.jl")
        print(result.stdout)
        for expected in (
            "testme_part_1.jl: all tests passed",
            "testme_part_2.jl: all tests passed",
            "testme_part_3.jl: all tests passed",
            "Advanced: checks passed",
            "Status: READY TO PACKAGE",
        ):
            if expected not in result.stdout:
                raise RuntimeError(f"check_submission.jl did not report: {expected}")
        written = {path.name for path in (student / "outputs").iterdir()}
        if written != REPORTS:
            raise RuntimeError(f"Unexpected reports in outputs/: {sorted(written ^ REPORTS)}")

        # Export only after validation succeeds, with a checksum beside the ZIP.
        output_directory.mkdir(parents=True, exist_ok=True)
        archive_path = output_directory / candidate.name
        shutil.copy2(candidate, archive_path)
        digest = hashlib.sha256(archive_path.read_bytes()).hexdigest()
        checksum_path = output_directory / f"{candidate.name}.sha256"
        checksum_path.write_text(f"{digest}  {candidate.name}\n", encoding="utf-8")

    print(f"Wrote tested solution archive: {archive_path}")
    return archive_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist")
    arguments = parser.parse_args()
    build_release(arguments.output_dir.resolve())
