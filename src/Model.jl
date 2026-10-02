"""
    reaction_index(model, id::AbstractString) -> Int

Find the reaction column identified by `id`.

# Arguments
- `model`: a named tuple whose `reactions` vector lists identifiers in matrix
  column order.
- `id::AbstractString`: the reaction identifier to locate.

# Returns and errors
Return its one-based column index as an `Int`. Throw `ArgumentError` if the
identifier is absent. Leave the model unchanged.
"""
function reaction_index(model, id::AbstractString)
    j = findfirst(==(id), model.reactions);
    isnothing(j) && throw(ArgumentError("unknown reaction: $id"));
    return j;
end

"""
    load_model(; glucose_limit=1.0) -> NamedTuple

Load the supplied deGFP model with the PS3 operating conditions.

# Arguments
- `glucose_limit::Real`: a finite, nonnegative glucose supply limit in mM/h;
  the default is 1.0. Read model files from the data directory beside src/.

# Returns
A named tuple containing newly loaded model arrays and fixed expression rates:
- `S`: stoichiometric matrix, with species in rows and reactions in columns.
- `reactions`, `metabolites`: identifiers in column and row order.
- `lower`, `upper`: reaction-flux bounds in mM/h.
- `amino_acids`, `amino_acid_uptake`: names and matching supply-column indices.
- `protein`, `translation`: column indices for protein output and translation
  initiation, respectively.
- `transcription_rate`, `translation_limit`: fixed expression rate and capacity
  in mM/h. Separate supply and removal reactions have nonnegative fluxes.

# Errors
Throw `ArgumentError` for an invalid glucose limit and `ErrorException` if the
matrix dimensions disagree with the identifiers. File-reading and parsing errors
propagate to the caller. The data files are read without being changed.
"""
function load_model(; glucose_limit::Real=1.0)
    isfinite(glucose_limit) && glucose_limit >= 0 ||
        throw(ArgumentError("glucose_limit must be finite and nonnegative"));
    # Read the matrix and labels in their stored row and column order -
    directory = joinpath(@__DIR__, "..", "data");
    S = readdlm(joinpath(directory, "stoichiometry.tsv"), '\t', Float64);
    rows = [split(line, '\t') for line in readlines(joinpath(directory, "reactions.tsv"))[2:end]];
    reactions = String[row[1] for row in rows];
    lower = [parse(Float64, row[2]) for row in rows];
    upper = [parse(Float64, row[3]) for row in rows];
    metabolites = readlines(joinpath(directory, "metabolites.txt"));
    size(S) == (length(metabolites), length(reactions)) || error("model dimensions do not agree");
    aa_rows = [split(line, '\t') for line in readlines(joinpath(directory, "amino_acids.tsv"))[2:end]];
    amino_acids = String[row[1] for row in aa_rows];
    lookup = Dict(id => j for (j, id) in enumerate(reactions));
    amino_acid_uptake = [lookup[row[2]] for row in aa_rows];

    # Apply Case 1: external amino-acid supply and internal synthesis are available -
    # Reaction names in the conditions table replace legacy numeric indices.
    for row in [split(line, '\t') for line in readlines(joinpath(directory, "conditions.tsv"))[2:end]]
        j = lookup[row[1]];
        lower[j] = parse(Float64, row[2]);
        upper[j] = parse(Float64, row[3]);
    end
    upper[lookup["M_glc_D_c_exchange"]] = Float64(glucose_limit);

    # Calculate fixed expression bounds from the original deGFP parameters -
    # These constants are evaluated before optimization; the optimization is linear.
    transcript_length = 683.0; # kinetic length parameter retained from the source
    # Transcription uses 25 nt/s, 75 nM polymerase, and 5 nM plasmid with a
    # 3.5 nM saturation constant. The factors 3600, 1e6, and 10/11 convert
    # seconds to hours, nM to mM, and apply the fixed promoter activity.
    transcription_rate = (25.0 / transcript_length) * (75.0 / 1e6) *
        (5.0 / (3.5 + 5.0)) * 3600.0 * (10.0 / 11.0); # mM/h
    mrna = transcription_rate / 5.2; # mM; the transcript degradation constant is 5.2/h
    # Translation uses polysome gain 10, 3 nt per amino acid, 2 amino acids/s,
    # 0.0016 mM ribosomes, and an mRNA saturation constant of 0.045 mM.
    translation_limit = (10.0 * 3.0 * 2.0 / transcript_length) * 3600.0 *
        0.0016 * mrna / (0.045 + mrna); # mM/h
    # Fix transcription and transcript degradation; cap translation initiation -
    for id in ("transcriptional_initiation_deGFP", "mRNA_degradation_deGFP")
        lower[lookup[id]] = transcription_rate;
        upper[lookup[id]] = transcription_rate;
    end
    translation = lookup["translation_initiation_deGFP"];
    upper[translation] = translation_limit;
    return (; S, reactions, metabolites, lower, upper, amino_acids,
        amino_acid_uptake, protein=lookup["PROTEIN_export_deGFP"], translation,
        transcription_rate, translation_limit);
end
