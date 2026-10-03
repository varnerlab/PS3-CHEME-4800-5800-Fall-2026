# PS3 grading rubric

PS3 is worth a maximum of **4 points** for the Standard track: Parts 1–2 and
questions 1–2 in [responses.md](responses.md). There are **48 Standard tests**,
24 each in [testme_part_1.jl](testme_part_1.jl) and
[testme_part_2.jl](testme_part_2.jl).

| Score | Requirements |
|:---:|---|
| 0 | Nothing runs; no tests pass. |
| 1 | Some Standard tests pass, but not a majority (1–24). |
| 2 | A majority of Standard tests pass, but at least one fails (25–47). |
| 3 | All 48 Standard tests pass, but questions 1–2 are not both satisfactorily answered. |
| 4 | All 48 Standard tests pass and questions 1–2 are satisfactorily answered. |

**Advanced track: one Magic Point.** Complete Part 3 and question 3 as well as
the Standard work. Earning the point requires a Standard score of 4, all 24
tests in [testme_part_3.jl](testme_part_3.jl) passing, the supply-range reports, and a
satisfactory answer to question 3. The teaching team reviews the code and
written responses before awarding the point, once for PS3. The Standard score
remains out of 4. Skipping or leaving Advanced work unfinished does not reduce it.

Select your track in `TRACK.txt` (`standard` or `advanced`) and run
`check_submission.jl`. On the Standard track, the checker ignores question 3's
placeholder and the `supply_ranges` starter. Passing checks does not
automatically award a Magic Point.

To grade, the instructor copies your `TRACK.txt`, `src/Compute.jl`, any helper
files, and `Include.jl` into a clean copy of the assignment, reruns the tests,
and then reads your `responses.md`.

The policies in [README.md](README.md) also apply: no initial submission, or
copying the reference solution, scores `0`.
Neither can earn a Magic Point. A revision of a qualifying initial submission
can still earn it after teaching-team review.
