# Calculation inputs

- [models/runs](models/runs): four all-locus configurations and three focused excerpts.
- [models/pedigrees](models/pedigrees): encoded parent-child relationships.
- [models/founder_scenarios](models/founder_scenarios): allele and explicit genotype priors.
- [models/identity_hypotheses](models/identity_hypotheses): identity-assignment input.
- [data](data): persons, reported observations, frequency data and locus-specific mapping.
- [Reported repeat counts](data/observed_str.reported_repeat_counts.json): reference transcription of the published STR calls; not laboratory raw data.
- [Observed-call coding](data/color_to_reported_str_allele.json): independently checked repeat-count interpretation of scored allele tokens.
- [Operational reference frequencies](data/allele_frequencies.reference.json): founder allele probabilities, source bins and zero-frequency floor.
- [Fixed background frequencies](data/background_allele_frequencies.json): shared by every run and kept fixed under founder variation.
- [Observation model](docs/OBSERVATION_MODEL_SPEC.md): normalisation and missing/partial calls.
- [Common evidence audit](docs/PRIMARY_EVIDENCE_AUDIT.md): scored sample-profile consistency.
- [Reference-frequency sources](docs/FREQUENCY_SOURCES.md): source URLs, selected populations and optional retrieval into an ignored cache.

Complete third-party webpages are not included. The calculations use the
prepared JSON; source retrieval is a separate provenance-audit step. See the
[licensing scope](../LICENSING.md) for third-party exclusions.

See the [main guide](../README.md#input-files-and-allele-coding) before changing
allele coding or interpreting a historical name in a JSON identifier.
