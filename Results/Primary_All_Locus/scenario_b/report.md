# Pedigree Likelihood Report

Run configuration: `Data/models/runs/scenario_b_primary_models_all_loci_generative_str.json`
Engine: `factors`
Observation model: `{"background_error": 0.001, "background_frequency_file": "../../data/background_allele_frequencies.json", "background_source": "fixed_reference_frequencies", "dropout": 0.2, "missing": "condition_on_call_available", "type": "generative_str"}`

## Ranking

| Rank | Model | logL | ΔlogL | Weight |
|---:|---|---:|---:|---:|
| 1 | hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity | -346.193537 | 0.000000 | 99.943% |
| 1 | belmonte_v1_scenario_b_all_loci / belmonte_v1_scenario_b_all_loci / amarna_no_identity | -346.193537 | 0.000000 | 99.943% |
| 2 | phizackerley_v2_scenario_b_all_loci / phizackerley_v2_scenario_b_all_loci / amarna_no_identity | -354.044595 | -7.851058 | 0.039% |
| 3 | variant_e1_scenario_b_all_loci / variant_e1_scenario_b_all_loci / amarna_no_identity | -355.104079 | -8.910542 | 0.013% |
| 4 | dodson_v2_scenario_b_all_loci / dodson_v2_scenario_b_all_loci / amarna_no_identity | -356.444069 | -10.250532 | 0.004% |
| 5 | tawfik_v1_scenario_b_all_loci / tawfik_v1_scenario_b_all_loci / amarna_no_identity | -358.203049 | -12.009512 | 0.001% |

## Support Ratios (LR for sensitivity settings; BF for population-prior comparisons)

| Comparison | Ratio |
|---|---:|
| hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity vs. belmonte_v1_scenario_b_all_loci / belmonte_v1_scenario_b_all_loci / amarna_no_identity | 1 |
| hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity vs. phizackerley_v2_scenario_b_all_loci / phizackerley_v2_scenario_b_all_loci / amarna_no_identity | 2568.45 |
| hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity vs. variant_e1_scenario_b_all_loci / variant_e1_scenario_b_all_loci / amarna_no_identity | 7409.68 |
| hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity vs. dodson_v2_scenario_b_all_loci / dodson_v2_scenario_b_all_loci / amarna_no_identity | 28297.6 |
| hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity vs. tawfik_v1_scenario_b_all_loci / tawfik_v1_scenario_b_all_loci / amarna_no_identity | 164310 |

## Locus Details

### Rank 1: hawass_v2_scenario_b_all_loci / hawass_v2_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -33.880439 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.230757 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.353854 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.926235 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.139146 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -59.617661 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.148378 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -28.897067 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 1: belmonte_v1_scenario_b_all_loci / belmonte_v1_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -33.880439 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.230757 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.353854 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.926235 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.139146 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -59.617661 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.148378 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -28.897067 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 2: phizackerley_v2_scenario_b_all_loci / phizackerley_v2_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -34.845952 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -65.514320 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.888273 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.395989 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.698349 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -61.697052 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.518123 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.486537 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 3: variant_e1_scenario_b_all_loci / variant_e1_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -34.845952 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -66.309692 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.888273 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.382043 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -36.401877 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -61.697052 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.092654 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.486537 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 4: dodson_v2_scenario_b_all_loci / dodson_v2_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.539008 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -66.528276 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -38.581412 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -49.382045 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -35.640271 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -61.591730 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.518577 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -29.662751 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |

### Rank 5: tawfik_v1_scenario_b_all_loci / tawfik_v1_scenario_b_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.539015 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -65.514291 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -38.581415 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -50.089127 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -36.391310 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -62.390193 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.518131 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -30.179568 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
