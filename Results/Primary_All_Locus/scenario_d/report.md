# Pedigree Likelihood Report

Run configuration: `Data/models/runs/scenario_d_primary_models_all_loci_generative_str.json`
Engine: `factors`
Observation model: `{"background_error": 0.001, "background_frequency_file": "../../data/background_allele_frequencies.json", "background_source": "fixed_reference_frequencies", "dropout": 0.2, "missing": "condition_on_call_available", "type": "generative_str"}`

## Ranking

| Rank | Model | logL | ΔlogL | Weight |
|---:|---|---:|---:|---:|
| 1 | hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity | -352.099526 | 0.000000 | 99.996% |
| 1 | belmonte_v1_scenario_d_all_loci / belmonte_v1_scenario_d_all_loci / amarna_no_identity | -352.099526 | 0.000000 | 99.996% |
| 2 | phizackerley_v2_scenario_d_all_loci / phizackerley_v2_scenario_d_all_loci / amarna_no_identity | -362.148757 | -10.049231 | 0.004% |
| 3 | tawfik_v1_scenario_d_all_loci / tawfik_v1_scenario_d_all_loci / amarna_no_identity | -366.307280 | -14.207754 | 0.000% |
| 4 | dodson_v2_scenario_d_all_loci / dodson_v2_scenario_d_all_loci / amarna_no_identity | -367.683344 | -15.583818 | 0.000% |
| 5 | variant_e1_scenario_d_all_loci / variant_e1_scenario_d_all_loci / amarna_no_identity | -374.026992 | -21.927466 | 0.000% |

## Support Ratios (LR for sensitivity settings; BF for population-prior comparisons)

| Comparison | Ratio |
|---|---:|
| hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity vs. belmonte_v1_scenario_d_all_loci / belmonte_v1_scenario_d_all_loci / amarna_no_identity | 1 |
| hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity vs. phizackerley_v2_scenario_d_all_loci / phizackerley_v2_scenario_d_all_loci / amarna_no_identity | 23138 |
| hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity vs. tawfik_v1_scenario_d_all_loci / tawfik_v1_scenario_d_all_loci / amarna_no_identity | 1.4803e+06 |
| hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity vs. dodson_v2_scenario_d_all_loci / dodson_v2_scenario_d_all_loci / amarna_no_identity | 5.86093e+06 |
| hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity vs. variant_e1_scenario_d_all_loci / variant_e1_scenario_d_all_loci / amarna_no_identity | 3.33409e+09 |

## Locus Details

### Rank 1: hawass_v2_scenario_d_all_loci / hawass_v2_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.303864 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus2 | -66.502441 | blue, cyan, green, other, red | 15 | max factor rows 15015 |
| Locus3 | -37.874963 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus4 | -50.378465 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus5 | -35.686264 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus6 | -58.255602 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus7 | -39.219495 | blue, green, other, red, yellow | 15 | max factor rows 15015 |
| Locus8 | -28.878432 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |

### Rank 1: belmonte_v1_scenario_d_all_loci / belmonte_v1_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.303864 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus2 | -66.502441 | blue, cyan, green, other, red | 15 | max factor rows 15015 |
| Locus3 | -37.874963 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus4 | -50.378465 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus5 | -35.686264 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus6 | -58.255602 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus7 | -39.219495 | blue, green, other, red, yellow | 15 | max factor rows 15015 |
| Locus8 | -28.878432 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |

### Rank 2: phizackerley_v2_scenario_d_all_loci / phizackerley_v2_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.935833 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -67.709028 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -39.641735 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -52.180848 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -36.139279 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -61.117339 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.869656 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -29.555039 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |

### Rank 3: tawfik_v1_scenario_d_all_loci / tawfik_v1_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -36.628908 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus2 | -67.709028 | blue, cyan, green, other, red | 15 | max factor rows 15015 |
| Locus3 | -40.334877 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus4 | -52.873990 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus5 | -36.832259 | blue, cyan, green, other, red, yellow | 21 | max factor rows 45291 |
| Locus6 | -61.810478 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |
| Locus7 | -39.869674 | blue, green, other, red, yellow | 15 | max factor rows 15015 |
| Locus8 | -30.248066 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 114982 |

### Rank 4: dodson_v2_scenario_d_all_loci / dodson_v2_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -36.628871 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -69.095306 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -40.334879 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -52.180850 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -36.696625 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -63.353178 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.869628 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -29.524007 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |

### Rank 5: variant_e1_scenario_d_all_loci / variant_e1_scenario_d_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.935833 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -78.706132 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -39.641735 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -52.180872 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -37.447076 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -61.117339 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.442965 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -29.555039 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
