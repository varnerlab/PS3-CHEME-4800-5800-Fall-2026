"""Build a student ZIP using an explicit file list; exclude instructor solutions."""
from pathlib import Path
import argparse
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = ["README.md", "RUBRIC.md", "responses.md", "Project.toml", "Manifest.toml",
              "Include.jl", "test_support.jl", "testme_part_1.jl", "testme_part_2.jl",
              "check_submission.jl", "runproduction.jl", "LICENSE", "THIRD_PARTY_NOTICES.txt"]
DIRECTORIES = ["src", "data"]

def student_files():
    files = [ROOT / name for name in ROOT_FILES]
    files += [p for name in DIRECTORIES for p in sorted((ROOT / name).rglob("*"))
              if p.is_file() and p.name != ".DS_Store" and "__pycache__" not in p.parts]
    assert all(p.is_file() for p in files), "A required student file is missing"
    return files

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--draft", action="store_true", help="Allow dates marked to be announced")
    args = parser.parse_args()
    if "To be announced" in (ROOT / "README.md").read_text() and not args.draft:
        parser.error("Set the release and submission dates before building a release, or use --draft")
    out = ROOT / "instructor" / "artifacts"
    out.mkdir(exist_ok=True)
    suffix = "-draft" if args.draft else ""
    archive = out / f"PS3-CHEME-4800-5800-Fall-2026-student{suffix}.zip"
    prefix = "PS3-CHEME-4800-5800-Fall-2026"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in student_files():
            bundle.write(path, str(Path(prefix) / path.relative_to(ROOT)))
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        assert not any(part in Path(name).parts for name in names
                       for part in ("solution", "instructor", ".git", "outputs"))
        compute = bundle.read(prefix + "/src/Compute.jl").decode()
        assert compute.count('error("Complete ') == 3, "Student starter must contain exactly three stubs"
    print(f"Built {archive} ({len(names)} files); no reference solution or generated answers included.")

if __name__ == "__main__":
    main()
