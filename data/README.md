# deGFP model data

Source: [Sequence-Specific-FBA-CFPS-Publication-Code](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code),
commit `e003d4b630f09f9ad0d177737cb13832a9071603`, deGFP model.

Model publication: Vilkhovoy, M.; Horvath, N.; Shih, C.-H.; Wayman, J. A.;
Calhoun, K.; Swartz, J.; Varner, J. D. (2018). **Sequence Specific Modeling of
E. coli Cell-Free Protein Synthesis.** *ACS Synthetic Biology*, **7**(8),
1844–1857. [DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465).
[PubMed: 29944340](https://pubmed.ncbi.nlm.nih.gov/29944340/).

| File | Contents |
|---|---|
| [stoichiometry.tsv](stoichiometry.tsv) | 146 × 265 stoichiometric matrix, without headers. |
| [reactions.tsv](reactions.tsv) | Reaction identifiers, original bounds, and equations in column order. |
| [metabolites.txt](metabolites.txt) | Species identifiers in row order. |
| [amino_acids.tsv](amino_acids.tsv) | Amino-acid names and their supply-reaction identifiers. |
| [conditions.tsv](conditions.tsv) | Named reaction-bound updates from the original Case 1. |
| [provenance.json](provenance.json) | Source commit and SHA-256 hashes of the original files. |

The matrix coefficients are unchanged. All 265 written reaction equations were
checked against the matrix during import. All original species-balance bounds
are zero. The source's nonnegative forward and reverse reaction columns are kept.

The loader applies the source's Case 1 (amino-acid uptake and synthesis), T7
promoter, and 5 nM plasmid setting. It preserves the original transcription and
translation formulas. PS3 changes glucose supply from 30 to **1 mM/h** and adds
the equal or shared amino-acid budget. Oxygen remains capped at **100 mM/h**.
The source's designated amino-acid degradation reactions and amino-acid removal
reactions are blocked. Every original individual amino-acid uptake ceiling is
30 mM/h before applying the student's allocation constraints.

The fixed transcription rate is approximately **0.005285 mM/h**. The fixed
translation limit is approximately **0.011176 mM/h**, or **11.176 µM/h**.
The expression parameters are supplied model inputs; PS3 does not ask students
to infer them from a protein sequence. The source uses a kinetic transcript-length
parameter of 683, while its written transcription reaction consumes 678
nucleotides. Both original quantities are retained for this teaching adaptation.

Protein output represents removal from the balanced network into the accumulating
product pool. The static optimization does not model substrate depletion,
changing enzyme activity, or final batch concentration. Rates use reaction-mixture
volume (mM/h), rather than a cell dry-weight normalization.

The original source files contain the MIT permission notice reproduced in
[THIRD_PARTY_NOTICES.txt](../THIRD_PARTY_NOTICES.txt).
