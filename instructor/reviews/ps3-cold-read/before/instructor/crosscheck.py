"""Compare all Julia curve values with an independent SciPy/HiGHS reconstruction.

Run: python3 instructor/crosscheck.py /path/to/pinned/publication/checkout
Requires NumPy and SciPy. These are instructor-only validation dependencies.
"""
from pathlib import Path
import csv
import hashlib
import json
import sys
import numpy as np
from scipy.optimize import linprog

root = Path(__file__).resolve().parents[1]
source = Path(sys.argv[1])
provenance = json.loads((root / "data/provenance.json").read_text())
for name, digest in provenance["source_sha256"].items():
    assert hashlib.sha256((source / name).read_bytes()).hexdigest() == digest
text = (source / "deGFP/DataDictionary.jl").read_text()
block = text.split("default_bounds_array = [")[1].split("];")[0]
bounds = np.array([[float(x) for x in line.split(";")[0].split()]
                   for line in block.splitlines() if ";" in line])
S = np.loadtxt(source / "deGFP/Network.dat")
# Reconstruct Bounds.jl Case 1 using the source's one-based reaction numbers.
upper_updates = {195: 100, 207: 1, 211: 0, 264: 0, 265: 0}
upper_updates.update({j: 0 for j in range(214, 218)})
upper_updates.update({j: 30 for j in range(224, 263, 2)})
upper_updates.update({j: 0 for j in range(225, 264, 2)})
upper_updates.update({j: 100 for j in range(76, 100)})
upper_updates.update({j: 0 for j in range(100, 117)})
for reaction, value in upper_updates.items():
    bounds[reaction - 1, 1] = value
tx = (25 / 683) * (75 / 1e6) * (5 / 8.5) * 3600 * (10 / 11)
mrna = tx / 5.2
tl = (10 * 3 * 2 / 683) * 3600 * .0016 * mrna / (.045 + mrna)
bounds[166] = tx
bounds[168] = tx
bounds[169, 1] = tl
c = np.zeros(265)
c[193] = -1
A = np.zeros((1, 265))
A[0, 223:262:2] = 1
rows = list(csv.DictReader((root / "instructor/reference/production-versus-budget.csv").open()))
errors = []
residuals = []
for row in rows:
    budget = float(row["budget_mM_per_h"])
    for strategy in ("equal", "optimized"):
        b = bounds.copy()
        kwargs = {}
        if strategy == "equal":
            b[223:262:2, 1] = np.minimum(b[223:262:2, 1], budget / 20)
        else:
            kwargs = {"A_ub": A, "b_ub": [budget]}
        result = linprog(c, A_eq=S, b_eq=np.zeros(146), bounds=b, method="highs", **kwargs)
        assert result.success, result.message
        errors.append(abs(-1000 * result.fun - float(row[strategy + "_uM_per_h"])))
        residuals.append(float(np.abs(S @ result.x).max()))
assert max(errors) < 1e-6, max(errors)
report = {"comparisons": len(errors), "maximum_rate_difference_uM_per_h": max(errors),
          "maximum_HiGHS_balance_residual_mM_per_h": max(residuals),
          "source_commit": provenance["commit"]}
directory = root / "instructor/validation-output"
directory.mkdir(exist_ok=True)
(directory / "independent-crosscheck.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
