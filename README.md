# Amarna Pedigree Analysis

Exact full-pedigree likelihood calculations for the published nuclear STR
profiles of the Amarna royal mummies. This is a self-contained computational
research archive: it contains Python code, documented JSON inputs, numerical
results and the checks needed to inspect and reproduce the analysis. No
article manuscript or journal PDF is included or required to use it.

The reported STR calls, operational founder frequencies and fixed background
distribution are documented separately under [Data](Data). The results use
the normalised dropout-and-background-error observation model described below.
The published calls are transcribed observations, not laboratory raw data.

The comparison evaluates **five unique scored genetic structures represented
by six labelled model forms**. Hawass-derived and Belmonte occupy the same
sampled inheritance positions and are normalised as **one** genetic structure.
Historical names are interpretations of those positions, not genetic outputs.

## Contents

- [Start here](#start-here)
- [Repository map](#repository-map)
- [What is being compared?](#what-is-being-compared)
- [Founder-background settings](#founder-background-settings)
- [Reading and interpreting the results](#reading-and-interpreting-the-results)
- [How the calculation works](#how-the-calculation-works)
- [Input files and allele coding](#input-files-and-allele-coding)
- [Reproducing the calculations](#reproducing-the-calculations)
- [Script reference](#script-reference)
- [Sources and retrieval](Data/docs/FREQUENCY_SOURCES.md)
- [Troubleshooting and verification](#troubleshooting-and-verification)
- [Citing this project](#citing-this-project)
- [Archiving and reuse](#archiving-and-reuse)

## Start here

**To read the findings without running Python:**

1. Read the [primary results summary](Results/Primary_All_Locus/summary/report.md#unique-scored-structure-rankings).
2. Open the [five-structure ranking CSV](Results/Primary_All_Locus/summary/structural_ranking_summary.csv).
3. Read the [all-locus robustness report](Results/Primary_All_Locus_Sensitivity/report.md).
4. Consult the [focused relationship stability report](Results/Focused_Robustness/stability_report.md).
5. Read the [observation-model specification](Data/docs/OBSERVATION_MODEL_SPEC.md) and [source retrieval guide](Data/docs/FREQUENCY_SOURCES.md).

**To verify the repository locally:** install Python 3.10 or later, clone the
repository, and run from its root:

```bash
git clone https://github.com/mintho/amarna-pedigree-analysis.git
cd amarna-pedigree-analysis
python3 tools/reproduce.py check
python3 tools/reproduce.py summary
```

Access to a private repository requires GitHub authentication. The scientific
scripts work with Python's standard library. Installing NumPy enables the
validated exact tensor-contraction acceleration used for the supplied results:

```bash
python3 -m pip install -r requirements-acceleration.txt
```

Without NumPy, the same factor engine uses the slower sparse implementation.
`check` runs unit tests, validates all seven configurations, checks input paths
and the common scored evidence, verifies distinct focused founder settings,
and checks archived structural weights. `summary` regenerates the ranking
tables from the archived likelihoods and compares the CSVs numerically.

## Repository map

| Location | Contents |
|---|---|
| [Code/src](Code/src) | Exact likelihood engines, genetics primitives, audits and report generators |
| [Code/tests](Code/tests) | Genetics and structural-normalisation tests |
| [Data/models/runs](Data/models/runs) | Four primary and three focused run configurations |
| [Data/models/pedigrees](Data/models/pedigrees) | Parent-child graph definitions |
| [Data/models/founder_scenarios](Data/models/founder_scenarios) | Locus frequencies and person-specific founder priors |
| [Data/data](Data/data) | Sample observations, person-node lists, allele mapping and frequency provenance |
| [Data/docs](Data/docs) | Observation-model specification and common-evidence audit |
| [Results/Primary_All_Locus](Results/Primary_All_Locus) | Scenario likelihoods and five-structure ranking tables |
| [Results/Primary_All_Locus_Sensitivity](Results/Primary_All_Locus_Sensitivity) | Full-ranking parameter grid and eight-locus jackknife |
| [Results/Focused_Robustness](Results/Focused_Robustness) | Three local pedigree robustness analyses |
| [Results/Frequency_prior_audit](Results/Frequency_prior_audit) | Reference-frequency extraction audit |
| [Results/Validation](Results/Validation) | Unit-test output, input audit and independent exact-engine comparisons |
| [Results/Manuscript_Tables/numerical_facts.json](Results/Manuscript_Tables/numerical_facts.json) | Machine-readable summary of the numerical findings |
| [Data/docs/FREQUENCY_SOURCES.md](Data/docs/FREQUENCY_SOURCES.md) | Source URLs, population choices and optional source retrieval |
| [tools/reproduce.py](tools/reproduce.py) | Reproduction commands with explicit analysis settings |
| [tools/verify_repository.py](tools/verify_repository.py) | Fast consistency and integrity checks |
| [tools/result_provenance.py](tools/result_provenance.py) | Code/input/result fingerprints for the reviewed calculation |

New calculations write to `reproduced/`, which is excluded from Git.
Archived inputs and results are not overwritten by the reproduction driver.

### Keeping inputs, code and results aligned

The canonical inputs are in `Data/`; the reviewed outputs are in `Results/`.
The [calculation provenance record](Results/Validation/calculation_provenance.json)
binds them to the calculation code. A changed JSON input or engine file makes
the integrity check fail until a new calculation has been reproduced and reviewed.

To review an intentional change, run `python3 tools/reproduce.py all --refresh`.
This writes separate outputs to `reproduced/` without accepting or promoting
them. Compare those outputs, archive the replaced local files outside this
repository, and promote only the reviewed outputs to `Results/`. The optional
`python3 tools/manuscript_numbers.py` export regenerates `numerical_facts.json`
and English/German LaTeX rows from those results. Its LaTeX output folders
are ignored by Git and are not required for numerical reproduction.

After fresh reproduction has confirmed the promoted outputs, explicitly record
them with `python3 tools/result_provenance.py --record-verified-calculation`.
Update the file inventory with `python3 tools/update_manifest.py`, then run
`check` and a normal numerical reproduction again. The optional source audit
also runs in `all` and needs archive access or an ignored local cache.
The provenance record proves file identity, not
scientific correctness; it must never be updated just to hide a discrepancy.

## What is being compared?

The primary comparison uses the same eight-locus observed sample profiles in
every candidate. It evaluates the entire connected pedigree, not a series of
independent pairwise tests.

| Model form in configuration/result IDs | Meaning |
|---|---|
| `hawass_v2_scenario_*` | Hawass-derived historical reading; KV55 interpreted as Akhenaten |
| `belmonte_v1_scenario_*` | Belmonte historical reading; KV55 interpreted as Smenkhkare |
| `dodson_v2_scenario_*` | Dodson-derived sampled structure |
| `phizackerley_v2_scenario_*` | Phizackerley-derived sampled structure |
| `tawfik_v1_scenario_*` | Tawfik sampled structure |
| `variant_e1_scenario_*` | Exploratory structural variant E1 |

The `v1`/`v2` strings are fixed input identifiers, not separate numerical
analyses to combine. The first two rows have identical genetic likelihoods:
their sampled profiles and inheritance positions coincide. Their historical
names differ, but they contribute only once to structural normalisation.

The common evidence includes Yuya, Thuya, Tiye, Amenhotep III, KV55, KV35YL,
Tutankhamun, KV21A and the two KV62 fetuses. **KV21B is not scored or placed by
the primary calculation.** Its published profile is retained in the reported
repeat-count transcription for source completeness, but not added to the
scored evidence. Historical source proposals are not extra ranked inputs.

See the [evidence audit](Data/docs/PRIMARY_EVIDENCE_AUDIT.md) for profile-level
checks. Mitochondrial and Y-chromosomal haplogroups are not inputs to these
nuclear-STR likelihoods.

## Founder-background settings

The four primary configurations are:
[a](Data/models/runs/scenario_a_primary_models_all_loci_generative_str.json),
[b](Data/models/runs/scenario_b_primary_models_all_loci_generative_str.json),
[c](Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json),
[d](Data/models/runs/scenario_d_primary_models_all_loci_generative_str.json).

| Setting | Interpretation | Reported support ratio |
|---|---|---|
| a | Maximally constrained Iuy-placeholder support; restrictive, data-informed sensitivity setting | Likelihood ratio (LR) |
| b | Partly constrained Iuy-placeholder support; restrictive, data-informed sensitivity setting | Likelihood ratio (LR) |
| c | Population-prior founder comparison, with Ay retained on the Yuya-Thuya background line | Bayes factor (BF) |
| d | Population-prior founder comparison, with Ay treated as an unrelated founder | Bayes factor (BF) |

The Iuy treatments in a/b are **allele-support/prior constraints**, not encoded
claims that she was related through one or both parents. Ay and Iuy are
historically motivated placeholder labels. The graph and explicit prior
tables, not a historical label alone, define the calculation.

All primary founders use the verified modern reference frequencies unless
the stated explicit Iuy-node sensitivity prior overrides them. These are
proxy assumptions, not allele frequencies estimated from the scored mummies.
The background-error distribution uses a separate, fixed reference file.

The focused robustness analyses use a different three-setting ladder:
baseline reference frequencies, a 50% blend towards uniform frequencies,
and restricted listed-allele support (`other` cap 0.01). These are three
distinct distributions, not the primary a-d scenarios and not literal
genealogical endogamy models. Their generation parameters are retained in the
[founder-variant manifests](Data/data).

## Reading and interpreting the results

### Primary ranking

The shared Hawass-derived/Belmonte structure is highest-likelihood in all four
settings. Its nearest non-tied alternative is Phizackerley-derived.

| Setting | Shared-structure weight | Advantage over nearest alternative | Ratio |
|---|---:|---:|---:|
| a | 99.914159% | 7.468524 log-likelihood units | LR 1,752.02 |
| b | 99.943460% | 7.851058 log-likelihood units | LR 2,568.45 |
| c | 99.978939% | 8.536285 log-likelihood units | BF 5,096.38 |
| d | 99.995594% | 10.049231 log-likelihood units | BF 23,137.98 |

For exact values use the [structural CSV](Results/Primary_All_Locus/summary/structural_ranking_summary.csv)
and [summary JSON](Results/Primary_All_Locus/summary/summary.json).

**Example, scenario c:** the shared structure has log likelihood
`-348.3634211380427`. Subtracting the nearest alternative's log likelihood gives
`8.536285134607908`; `exp(8.536285134607908)` is approximately `5096.38`. Under
these assumptions, the published observations are about 5,096 times more
probable under the shared structure than under that alternative. This is not
a 5,096-fold probability that KV55's historical name is correct.

### Historical aliases and structural weights

The raw [scenario c JSON](Results/Primary_All_Locus/scenario_c/results.json)
and [label-level CSV](Results/Primary_All_Locus/summary/ranking_summary.csv)
retain six labelled outputs, but normalisation uses five unique `structure_id`
values. The two historical readings share the same rank and structural weight.
They are aliases of one scored structure: do not add their repeated weights
or treat them as independent evidence.

The authoritative [structural CSV](Results/Primary_All_Locus/summary/structural_ranking_summary.csv)
collapses the tied readings and normalises over five unique structures with
equal structural prior weights:

```text
delta_i = logL_i - max(logL)
weight_i = 100 * exp(delta_i) / sum_j exp(delta_j)
support(best, alternative) = exp(logL_best - logL_alternative)
```

A structural weight is conditional on this finite candidate set, the equal
structure priors, the reported observations and the modelling assumptions.
It is not an unconditional probability that a reconstruction is true, nor
evidence authenticating the DNA. Adding different candidate structures or
changing their priors can change the normalised weights.

### JSON and CSV fields

| Field | How to read it |
|---|---|
| `pedigree`, `founder_scenario`, `identity_hypothesis` | Input IDs identifying graph, priors and identity treatment |
| `log_likelihood` | Natural logarithm of the probability of the observed data under that model; less negative is better |
| `delta_log_likelihood` in ranking files | Candidate minus best; zero for best, negative for alternatives |
| `structure_id`, `weight_scope` | Identifies the scored structure and makes its normalisation explicit |
| `rank` | Structural rank; historical aliases share a rank |
| `weight_percent` | Normalisation over unique structures, repeated for historical aliases |
| `loci` | Locus-specific likelihoods and, for the factor engine, computational diagnostics |
| `observation_model` | Dropout/error settings applied to that calculation |
| `nearest_non_tied_bayes_factor` | Stored support-ratio field; interpret as LR in a/b and BF in c/d |
| `run_config`, `created_at_utc` | Provenance metadata, not quantities expected to match byte-for-byte in a rerun |

Support-ratio fields use the scenario-specific LR/BF interpretation above.

### Per-locus contributions

The [per-locus CSV](Results/Primary_All_Locus/summary/per_locus_vs_nearest.csv)
decomposes the shared structure's total advantage over its nearest alternative.
Here `delta_log_likelihood` is **best minus alternative** at each locus, so
positive values favour the shared structure; this sign convention differs
from the ranking CSV. Sum the eight contributions to recover the overall gap.
This is not eight independent historical-model rankings.

### Robustness outputs

| Design | What varies | Number of calculations |
|---|---|---:|
| [Full-ranking parameter grid](Results/Primary_All_Locus_Sensitivity/parameter_grid.csv) | Four scenarios x three error rates x three dropout rates | 36 |
| [Full-ranking leave-one-locus jackknife](Results/Primary_All_Locus_Sensitivity/leave_one_locus.csv) | Each of eight loci omitted under each scenario | 32 |
| [Focused matrices](Results/Focused_Robustness) | Three distinct focused founder priors x three error rates x three dropout rates | 27 per excerpt |
| [Focused jackknife](Results/Focused_Robustness/jackknife_report.md) | Omit L1 or L2 in each two-locus excerpt | 2 per excerpt |

The shared structure remains top-ranked in every archived full-ranking
perturbation and omission. All three focused excerpts retain their baseline
winner at 100% of grid points (0/27 rank flips) and in both locus omissions.
A rank flip means
that another candidate becomes highest-likelihood at that tested setting.
**These percentages describe the chosen grid, not posterior probabilities.**

The rankings are stable but their weight changes differ. Omitting L2 gives
relationship weights of 59.897% (Yuya-Thuya-Tiye), 97.960% (the immediate
parental excerpt) and 91.183% (Yuya-KV55). These compare a relationship with
an independent-person alternative, not with all five historical structures.
See the [focused jackknife report](Results/Focused_Robustness/jackknife_report.md)
and [numerical facts](Results/Manuscript_Tables/numerical_facts.json).
Neither test identifies KV55 as Akhenaten/Smenkhkare or KV35YL as
Baketaten: names that do not change sampled positions do not change likelihoods.

## How the calculation works

1. Load persons, reported observations, pedigree edges and founder priors from
   the run's referenced JSON files.
2. At each locus, construct possible unordered diploid genotypes. Form founder
   genotype probabilities from locus frequencies (Hardy-Weinberg proportions)
   unless an explicit person-specific prior overrides them.
3. Apply Mendelian transmission along the specified parent-child edges.
4. Apply the normalised `generative_str` observation probability to each
   available reported call.
5. Sum over unobserved genotype assignments. The factor engine does this by
   exact variable elimination; the reference engine enumerates assignments.
6. Multiply locus likelihoods, implemented as addition of natural logarithms,
   then rank candidates and form support ratios.
7. Collapse the two identical scored readings before structural normalisation.

No random sampling is used by these calculations, so there is no random seed
to set. Floating-point roundoff and metadata can differ between environments.
The driver checks numerical tables using absolute tolerance `1e-8` and relative
tolerance `1e-9`.

### Observation model

For the primary comparison, per-copy dropout `d = 0.2` and background-error
mixture `eta = 0.001` are modelling assumptions, not measured sample-specific
estimates. Dropout acts on two physical copies, including homozygous genotypes.
Authentic call probabilities are conditioned on at least one copy surviving.

```text
P(call | genotype) = (1 - eta) * P_authentic(call | genotype, d)
                    + eta * P_background(call | d)
P(missing call | genotype) = 1
```

The background-call distribution is fixed by the baseline frequency treatment,
not by a candidate's person-specific founder restrictions. Error-rate and
dropout perturbations deliberately change observation assumptions. Read the
[full specification](Data/docs/OBSERVATION_MODEL_SPEC.md) for the call state
space, normalisation, homozygotes and the background mixture.

Missing allele calls are not completed by local Mendelian deductions.
Incomplete profiles remain incomplete and their latent genotypes are summed over.

## Input files and allele coding

Run configurations refer to their inputs **relative to the configuration file**.
Keep `Data/data` and `Data/models` together. Paths in frequency-source metadata
and founder-variant manifests are relative to `Data/`; they describe provenance
rather than additional run configurations.

| JSON component | Meaning |
|---|---|
| `persons.*.json` | Person-node identifiers available to the pedigree |
| `observed_str.*.json` | Reported per-person, per-locus calls |
| `models/pedigrees/*.json` | Explicit parent-child relationships |
| `models/founder_scenarios/*.json` | Locus allele frequencies and optional explicit founder genotype priors |
| `models/identity_hypotheses/*.json` | Identity treatment referenced by the likelihood engine |
| `models/runs/*.json` | Links the inputs and selects the observation model |

A call with `alleles: []` is missing, one allele is partial, and two alleles
are a complete unordered call. Do not treat a partial single-allele call as
a reported homozygote or replace missing alleles with diagram implications.

The code uses `Locus1`-`Locus8`; their short labels are L1-L8:

| Internal ID | Short label | STR marker |
|---|---|---|
| `Locus1` | L1 | D13S317 |
| `Locus2` | L2 | D7S820 |
| `Locus3` | L3 | D2S1338 |
| `Locus4` | L4 | D21S11 |
| `Locus5` | L5 | D16S539 |
| `Locus6` | L6 | D18S51 |
| `Locus7` | L7 | CSF1PO |
| `Locus8` | L8 | FGA |

Allele tokens such as `blue`, `green` and `cyan` are **locus-specific states**.
The same token does not denote the same repeat count at different markers.
Use the [mapped allele key](Data/data/color_to_str_allele.duesseldorf_egypt.json)
and [documented mapping choices](Data/data/color_to_str_allele.manual_overrides.json),
not colour names alone, to recover repeat counts.

The [reported repeat counts](Data/data/observed_str.reported_repeat_counts.json)
are a transcription of Hawass et al. 2010, figure 1, p. 641, not laboratory raw
data. The [observed-call key](Data/data/color_to_reported_str_allele.json)
independently verifies the operational coding against all ten scored profiles.
The [operational frequency input](Data/data/allele_frequencies.reference.json)
supplies founder priors. It distinguishes actual reported repeat counts from
grouped source-table frequency bins, for example `6/6.2/6.3` at D7S820.
The [fixed background frequencies](Data/data/background_allele_frequencies.json)
are never changed by a founder-sensitivity variant.

The [extracted reference frequencies](Data/data/allele_frequencies.duesseldorf_egypt.json)
retain the population/source details: Egyptian frequencies where available,
with the documented Israeli fallback at D2S1338. The
[source guide](Data/docs/FREQUENCY_SOURCES.md) links each archived source page
and explains optional retrieval into an ignored local cache. Complete HTML
pages are not distributed. Displayed zeros receive
the stated floor `1e-5` in founder-prior construction; `other` represents
remaining probability mass as a latent category, not a named observed allele.

Historical names in person IDs identify positions in a represented reading.
They are not independently authenticated mummy identifications.

## Reproducing the calculations

Use these commands from the repository root. All driver tasks write generated
files beneath `reproduced/` and leave `Data/` and `Results/` intact. Full
pedigree calculations can be substantially slower and more memory-intensive
than validation or the two-locus excerpts; inspect the archived results first
if you only need to read them. Do not use brute-force enumeration for the full
eight-locus model set.

### 1. Validate inputs and tests

```bash
python3 tools/reproduce.py check
```

This checks the current archived material without recalculating all likelihoods.
It also writes a fresh common-evidence audit to `reproduced/validation/`.

### 2. Rebuild tables from stored likelihoods

```bash
python3 tools/reproduce.py summary
```

This is a fast table-generation check, not a fresh genetic calculation. It
regenerates and verifies the structural, label-level and per-locus CSVs in
`reproduced/Primary_All_Locus/summary/`.

### 3. Recalculate the four primary settings

```bash
python3 tools/reproduce.py primary
```

Runs all six labelled model forms at all eight loci under each of a-d using
the factor engine, then rebuilds the five-structure tables and compares them
against archived CSVs. Each `reproduced/Primary_All_Locus/scenario_*/` contains
`results.json` and `report.md`.

To run just scenario c directly:

```bash
python3 Code/src/compare.py \
  Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json \
  --engine factors --out_dir reproduced/Primary_All_Locus/scenario_c
```

### 4. Recalculate the focused robustness analyses

```bash
python3 tools/reproduce.py focused
```

Runs all three excerpts with **three distinct founder priors**, error rates
`0.0001, 0.001, 0.01` and dropout values `0, 0.1, 0.2` (27 points each).
The focused baseline uses dropout **0**, unlike the primary baseline **0.2**.
It then omits L1 and L2 separately and generates stability, jackknife and
classification reports. Numerical report tables are compared against the
archive. Per-excerpt `matrix.json` and `jackknife.json` are written alongside
the generated reports.

### 5. Recalculate the full-ranking grid and jackknife

```bash
python3 tools/reproduce.py sensitivity
```

This performs 36 complete comparison grids plus 32 leave-one-locus comparisons:
primary dropout values `0.1, 0.2, 0.3`, error rates `0.0001, 0.001, 0.01`,
then eight omissions under each scenario at baseline dropout/error. Output is
`reproduced/Primary_All_Locus_Sensitivity/`; both CSVs are checked against the
archive. This is the most extensive reproduction task.

### 6. Re-extract frequency provenance

```bash
python3 tools/reproduce.py priors
```

Downloads the archived reference pages on first use, caches them under
`reproduced/source_cache/duesseldorf_archive/`, and uses the recorded mapping
choices to produce an audit and allele-key JSON under
`reproduced/Frequency_prior_audit/`. Both directories are ignored by Git.
It does not replace operational priors. This optional source audit depends
on archive availability; numerical calculations use the included JSON and
do not require the webpages. See the [retrieval guide](Data/docs/FREQUENCY_SOURCES.md).

### 7. Run every task

```bash
python3 tools/reproduce.py all
```

Runs validation, frequency extraction, primary calculations, focused analyses
and full-ranking sensitivity in that order. Outputs from repeated invocations
replace files in `reproduced/`, never the archived results.
First-time source retrieval requires network access. To reproduce only the
numerical analysis offline, run `check`, `summary`, `primary`, `focused` and
`sensitivity` separately.

## Script reference

| Script/module | Responsibility |
|---|---|
| [genetics.py](Code/src/genetics.py) | Genotype states, founder probabilities, gametes, Mendelian transmission, normalised dropout/background call probabilities |
| [validate.py](Code/src/validate.py) | Validate run/input schemas, references and model consistency |
| [likelihood_factors.py](Code/src/likelihood_factors.py) | Exact sparse factor construction and variable elimination; production engine |
| [likelihood_enumeration.py](Code/src/likelihood_enumeration.py) | Exact assignment enumeration; reference engine for tractable examples; shared input handling |
| [compare.py](Code/src/compare.py) | Evaluate a run, rank labelled outputs and write likelihood reports |
| [summarize_primary_results.py](Code/src/summarize_primary_results.py) | Merge equivalent readings, normalise unique structures, generate ranking and per-locus tables |
| [primary_generative_sensitivity.py](Code/src/primary_generative_sensitivity.py) | The documented 36 full-ranking parameter points and 32 locus omissions |
| [sensitivity_matrix.py](Code/src/sensitivity_matrix.py) | Focused grid over founder priors, dropout and background error; use explicit grids via the driver |
| [locus_jackknife.py](Code/src/locus_jackknife.py) | Evaluate a baseline and leave each available locus out in turn |
| [stability_report.py](Code/src/stability_report.py) | Summarise grid winners, rank flips and changes in weights |
| [locus_jackknife_report.py](Code/src/locus_jackknife_report.py) | Summarise locus-omission results |
| [robustness_conclusions.py](Code/src/robustness_conclusions.py) | Classify the focused sensitivity patterns and write interpretive reports |
| [background_error_sensitivity.py](Code/src/background_error_sensitivity.py) | Single-axis error-rate helper used by the sensitivity toolchain |
| [founder_prior_sensitivity.py](Code/src/founder_prior_sensitivity.py) | Founder-prior helper used by the sensitivity toolchain |
| [benchmark.py](Code/src/benchmark.py) | Time and compare enumeration/factor likelihood engines on tractable configurations |
| [extract_duesseldorf_egypt_priors.py](Code/src/extract_duesseldorf_egypt_priors.py) | Retrieve and parse archived source tables; record frequencies/provenance |
| [build_color_str_mapping.py](Code/src/build_color_str_mapping.py) | Map locus-specific visual tokens to repeat counts, with recorded overrides |
| [build_duesseldorf_founder_scenario.py](Code/src/build_duesseldorf_founder_scenario.py) | Form mapped reference-frequency founder inputs with the documented zero floor |
| [define_founder_prior_variants.py](Code/src/define_founder_prior_variants.py) | Prepare broad/restricted focused founder-prior variants and manifests |
| [audit_primary_evidence.py](Code/src/audit_primary_evidence.py) | Check common scored profiles, fetus-parent assignments and exclusion of KV21B |
| [prepare_inputs.py](tools/prepare_inputs.py) | Synchronise operational primary founder frequencies and fixed background with the verified reference key |
| [validate_engines.py](tools/validate_engines.py) | Compare exact sparse, accelerated and enumeration implementations on current inputs |
| [manuscript_numbers.py](tools/manuscript_numbers.py) | Optional numerical-facts and bilingual LaTeX export from Results |
| [update_manifest.py](tools/update_manifest.py) | Inventory the portable bundle, excluding retrieved pages and local publication assets |

Input preparation helpers are included for provenance and inspection; use the
recorded inputs to reproduce this analysis. Their generic defaults need not match
the retained grids or repository-relative paths. The reproduction driver passes
the actual paths and numerical settings explicitly.

## Troubleshooting and verification

| Symptom | Check |
|---|---|
| `FileNotFoundError` for a configuration/input | Run from the repository root and retain the directory layout; use the driver |
| `ModuleNotFoundError: src` when running tests directly | Use `tools/reproduce.py check`, which sets the test import path |
| Two entries with the same top weight | Historical aliases repeat one structural weight; use the five-structure CSV |
| Ratios called BF in a raw sensitivity report | Interpret a/b as descriptive LRs; c/d are the formal BF comparisons |
| Different timestamps or local paths after rerunning | Expected provenance differences; compare numerical fields, not whole-file bytes |
| Slow calculation/high memory use | Install NumPy for exact contraction acceleration; avoid full-pedigree enumeration |
| Different grid percentages | Confirm exactly three distinct focused priors and the explicit error/dropout grids |
| A single coloured token seems inconsistent across markers | Mapping is locus-specific; consult the mapped allele key |

For a small independent-engine comparison:

```bash
python3 Code/src/benchmark.py \
  Data/models/runs/focused_yuya_thuya_tiye_trio_l1_l2_run.json \
  --engine both --repeats 1
```

See [test_genetics.py](Code/tests/test_genetics.py) for coverage of inheritance,
founder distributions, observation normalisation, complete/partial calls,
homozygous dropout and background errors. The additional
[repository tests](Code/tests/test_repository.py) check five-structure
normalisation and the reported support ratio.

`MANIFEST.sha256` records the distributed scientific inputs, scripts, results,
and documentation. `check` verifies these hashes when the manifest
is present. A new computation is expected to match numerically, not to have
the same provenance metadata or byte hash as an archived output.

Result provenance uses repository-relative paths. The seven executable
configurations are under `Data/models/runs/`.

## Citing this project

Author: Thomas Minzenmay, [ORCID 0009-0007-1730-4528](https://orcid.org/0009-0007-1730-4528).
The [CITATION.cff](CITATION.cff) file supplies GitHub's citation metadata.
The matching [.zenodo.json](.zenodo.json) supplies the deposit metadata used
by the Zenodo integration; Zenodo gives that file precedence over the CFF.
Version `1.0.0` is archived on Zenodo:
[10.5281/zenodo.23064916](https://doi.org/10.5281/zenodo.23064916).

Cite: Minzenmay, Thomas (2026). *Amarna Pedigree Analysis* (version 1.0.0)
[Software]. Zenodo. https://doi.org/10.5281/zenodo.23064916

Use this version-specific DOI for the archived calculations rather than the
mutable `main` branch. The [concept DOI](https://doi.org/10.5281/zenodo.23064915)
represents the project across versions. A later software version should receive
its own release metadata and version-specific citation.

## Archiving and reuse

This computational archive can be inspected and reproduced independently of
an article. The GitHub repository is connected to Zenodo. The published
archive was verified file-for-file against release `v1.0.0`; the release tag
and deposited files remain unchanged when citation documentation is updated.

The Python software is licensed under [MIT](LICENSE). Original documentation
and protectable contributions to result presentations are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), with attribution to
Thomas Minzenmay. Third-party observations and reference-frequency data are
explicitly excluded from these grants; see [licensing scope](LICENSING.md).
Complete source webpages, article PDFs and publication artwork are not
distributed in this repository. Source attribution and retrieval instructions
are retained alongside the extracted values.
