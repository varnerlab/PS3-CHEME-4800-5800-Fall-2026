"""Compare all Julia curve values and supply ranges with an independent SciPy/HiGHS reconstruction.

Run: python3 instructor/crosscheck.py /path/to/pinned/publication/checkout
Requires NumPy and SciPy. These are instructor-only validation dependencies.
Read the pinned source files and the reference production, supply-range, and
range-width tables in instructor/reference/. Write or overwrite
instructor/validation-output/independent-crosscheck.json and print the comparison
summary. Hash mismatches, unsuccessful solves, rate differences of 1e-6 μM/h or
more, range-endpoint or total-width differences of 1e-6 mM/h or more, or different
fixed-supply counts raise AssertionError. Source data are unchanged.
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
    bounds[reaction - 1, 1] = value  # convert source reaction numbers to Python indices
# Reproduce the fixed transcription rate and translation capacity from the source.
tx = (25 / 683) * (75 / 1e6) * (5 / 8.5) * 3600 * (10 / 11)
mrna = tx / 5.2
tl = (10 * 3 * 2 / 683) * 3600 * .0016 * mrna / (.045 + mrna)
bounds[166] = tx
bounds[168] = tx
bounds[169, 1] = tl
c = np.zeros(265)
c[193] = -1  # SciPy minimizes; a negative coefficient maximizes protein output
A = np.zeros((1, 265))
A[0, 223:262:2] = 1
rows = list(csv.DictReader((root / "instructor/reference/production-versus-budget.csv").open()))
supply = list(range(223, 262, 2))  # zero-based supply columns, in amino-acid table order

def scenario(budget, strategy):
    """Return the (bounds, inequality kwargs) pair for one budget and strategy.

    budget is in mM/h and strategy is "equal" or "optimized". Copy the Case 1
    bounds; leave the module arrays unchanged.
    """
    b = bounds.copy()
    if strategy == "equal":
        b[supply, 1] = np.minimum(b[supply, 1], budget / 20)
        return b, {}
    return b, {"A_ub": A, "b_ub": [budget]}

def ranges(budget, strategy, fraction):
    """Return a list of (minimum, maximum) supply fluxes in mM/h by flux variability analysis.

    Maximize protein output, require at least fraction times that optimum, then
    minimize and maximize each supply flux with HiGHS. Raise AssertionError for
    an unsuccessful solve. Inputs and module arrays are unchanged.
    """
    b, kwargs = scenario(budget, strategy)
    first = linprog(c, A_eq=S, b_eq=np.zeros(146), bounds=b, method="highs", **kwargs)
    assert first.success, first.message
    z = -first.fun  # maximum protein output in mM/h
    # c holds -1 at the protein column, so the row c @ v <= -fraction*z requires
    # protein output of at least fraction*z, with no tolerance, as in the assignment.
    A_fva = np.vstack([kwargs.get("A_ub", np.zeros((0, 265))), c.reshape(1, -1)])
    b_fva = list(kwargs.get("b_ub", [])) + [-fraction * z]
    out = []
    for j in supply:
        e = np.zeros(265)
        e[j] = 1.0
        lo = linprog(e, A_ub=A_fva, b_ub=b_fva, A_eq=S, b_eq=np.zeros(146), bounds=b, method="highs")
        hi = linprog(-e, A_ub=A_fva, b_ub=b_fva, A_eq=S, b_eq=np.zeros(146), bounds=b, method="highs")
        assert lo.success and hi.success, (budget, strategy, fraction, lo.message, hi.message)
        out.append((lo.fun, -hi.fun))
    return out

errors = []
residuals = []
# Compare production rates, allowing different optimal internal flux vectors.
for row in rows:
    budget = float(row["budget_mM_per_h"])
    for strategy in ("equal", "optimized"):
        b, kwargs = scenario(budget, strategy)
        result = linprog(c, A_eq=S, b_eq=np.zeros(146), bounds=b, method="highs", **kwargs)
        assert result.success, result.message
        errors.append(abs(-1000 * result.fun - float(row[strategy + "_uM_per_h"])))
        residuals.append(float(np.abs(S @ result.x).max()))
assert max(errors) < 1e-6, max(errors)

# Compare the unique supply-range endpoints at U=1 with the reference table.
table = list(csv.DictReader((root / "instructor/reference/supply-ranges.csv").open()))
range_errors = []
for label, strategy, fraction in (("equal", "equal", 1.0), ("optimized", "optimized", 1.0),
                                  ("optimized_99pct", "optimized", 0.99)):
    for (lo, hi), entry in zip(ranges(1.0, strategy, fraction), table, strict=True):
        range_errors.append(abs(lo - float(entry[label + "_min_mM_per_h"])))
        range_errors.append(abs(hi - float(entry[label + "_max_mM_per_h"])))
assert max(range_errors) < 1e-6, max(range_errors)

# Compare total range widths and fixed-supply counts at every plotted budget.
width_errors = []
for entry in csv.DictReader((root / "instructor/reference/range-width-versus-budget.csv").open()):
    budget = float(entry["budget_mM_per_h"])
    for strategy in ("equal", "optimized"):
        widths = [hi - lo for lo, hi in ranges(budget, strategy, 1.0)]
        width_errors.append(abs(sum(widths) - float(entry[strategy + "_total_width_mM_per_h"])))
        fixed = sum(w <= 1e-5 for w in widths)
        assert fixed == int(entry[strategy + "_fixed_supplies"]), (budget, strategy, fixed)
assert max(width_errors) < 1e-6, max(width_errors)

report = {"comparisons": len(errors), "maximum_rate_difference_uM_per_h": max(errors),
          "maximum_HiGHS_balance_residual_mM_per_h": max(residuals),
          "range_endpoint_comparisons": len(range_errors),
          "maximum_range_endpoint_difference_mM_per_h": max(range_errors),
          "width_comparisons": len(width_errors),
          "maximum_total_width_difference_mM_per_h": max(width_errors),
          "source_commit": provenance["commit"]}
directory = root / "instructor/validation-output"
directory.mkdir(exist_ok=True)
(directory / "independent-crosscheck.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
