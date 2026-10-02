"""Extract the pinned publication model into transparent PS3 data tables.

Run: python3 instructor/import_model.py /path/to/publication/checkout
Only the Python standard library is required. Read the source checkout without
changing it; overwrite the generated files in data/. Verify the pinned revision
and the agreement between written reactions and matrix coefficients first.
Failed checks raise AssertionError; file and parsing errors propagate.
"""
from pathlib import Path
import hashlib
import json
import re
import subprocess
import sys

COMMIT = "e003d4b630f09f9ad0d177737cb13832a9071603"
root = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1]).resolve()
revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
assert revision == COMMIT, f"Expected {COMMIT}, found {revision}"
text = (source / "deGFP/DataDictionary.jl").read_text()

def block(name):
    """Return the contents of the named Julia array assignment in the source text.

    name is an assignment identifier in the loaded deGFP DataDictionary.jl.
    Strip the opening 'name = [' and closing '];' markers. Leave the global
    source text unchanged; a missing opening marker raises IndexError.
    """
    return text.split(name + " = [", 1)[1].split("];", 1)[0]

reactions = re.findall(r'"([^"\n]+)"', block("list_of_reaction_strings"))
metabolites = re.findall(r'"([^"\n]+)"', block("list_of_metabolite_symbols"))
bounds = [[float(x) for x in line.split(";")[0].split()]
          for line in block("default_bounds_array").splitlines() if ";" in line]
species_bounds = [[float(x) for x in line.split(";")[0].split()]
                  for line in block("species_bounds_array").splitlines() if ";" in line]
matrix = [[float(x) for x in line.split()] for line in (source / "deGFP/Network.dat").read_text().splitlines()]
assert len(matrix) == len(metabolites) == 146
assert len(reactions) == len(bounds) == 265
assert all(len(row) == 265 for row in matrix)
assert all(pair == [0.0, 0.0] for pair in species_bounds)

# Independently reconstruct every coefficient from the written reaction equations.
reconstructed = [[0.0] * len(reactions) for _ in metabolites]
indices = {name: i for i, name in enumerate(metabolites)}
for j, reaction in enumerate(reactions):
    _, equation = reaction.split("::")
    left, right = equation.split(" --> ")
    for sign, side in [(-1, left), (1, right)]:
        for term in side.split("+"):
            if term == "[]":
                continue
            parts = term.split("*")
            coefficient = float(parts[0]) if len(parts) == 2 else 1.0
            reconstructed[indices[parts[-1]]][j] += sign * coefficient
assert reconstructed == matrix, "Reaction equations and matrix differ"
data = root / "data"
data.mkdir(exist_ok=True)
(data / "stoichiometry.tsv").write_text("".join("\t".join(f"{x:g}" for x in row) + "\n" for row in matrix))
(data / "metabolites.txt").write_text("\n".join(metabolites) + "\n")
ids = [reaction.split("::")[0] for reaction in reactions]
assert len(set(ids)) == len(ids)
lines = ["reaction\tlower\tupper\tequation\n"]
for reaction, (lo, hi) in zip(reactions, bounds):
    name, equation = reaction.split("::")
    lines.append(f"{name}\t{lo:g}\t{hi:g}\t{equation}\n")
(data / "reactions.tsv").write_text("".join(lines))
aa_names = ["Alanine", "Arginine", "Asparagine", "Aspartate", "Cysteine", "Glutamate",
            "Glutamine", "Glycine", "Histidine", "Isoleucine", "Leucine", "Lysine",
            "Methionine", "Phenylalanine", "Proline", "Serine", "Threonine",
            "Tryptophan", "Tyrosine", "Valine"]
(data / "amino_acids.tsv").write_text("amino_acid\tuptake_reaction\n" + "".join(
    f"{name}\t{ids[j]}\n" for name, j in zip(aa_names, range(223, 262, 2))))

# Preserve Case 1 bounds, including secretion restrictions. These are zero-based
# Python column indices; the original Julia reaction numbers are one larger.
conditions = {194: 100.0, 206: 30.0, 210: 0.0, 263: 0.0, 264: 0.0}
conditions.update({j: 0.0 for j in range(213, 217)})
conditions.update({j: 30.0 for j in range(223, 262, 2)})
conditions.update({j: 0.0 for j in range(224, 263, 2)})
conditions.update({j: 100.0 for j in range(75, 99)})
conditions.update({j: 0.0 for j in range(99, 116)})
(data / "conditions.tsv").write_text("reaction\tlower\tupper\n" + "".join(
    f"{ids[j]}\t0\t{value:g}\n" for j, value in sorted(conditions.items())))
files = ["deGFP/Network.dat", "deGFP/DataDictionary.jl", "deGFP/TXTLDictionary.jl", "Model/Bounds.jl", "Model/Solve.jl"]
provenance = {"repository": "https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code",
              "publication": {
                  "title": "Sequence Specific Modeling of E. coli Cell-Free Protein Synthesis",
                  "authors": ["Michael Vilkhovoy", "Nicholas Horvath", "Che-Hsiao Shih",
                              "Joseph A. Wayman", "Kara Calhoun", "James Swartz", "Jeffrey D. Varner"],
                  "year": 2018, "journal": "ACS Synthetic Biology", "volume": 7,
                  "issue": 8, "pages": "1844-1857", "doi": "10.1021/acssynbio.7b00465",
                  "pmid": "29944340", "url": "https://pubmed.ncbi.nlm.nih.gov/29944340/"},
              "commit": COMMIT, "source_sha256": {
                  f: hashlib.sha256((source / f).read_bytes()).hexdigest() for f in files}}
(data / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print("Imported 146 species and 265 reactions; all written equations match the matrix.")
