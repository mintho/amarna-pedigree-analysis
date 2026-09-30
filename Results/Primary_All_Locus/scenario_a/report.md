# Pedigree Likelihood Report

Run configuration: `Data/models/runs/scenario_a_primary_models_all_loci_generative_str.json`
Engine: `factors`
Observation model: `{"background_error": 0.001, "background_frequency_file": "../../data/background_allele_frequencies.json", "background_source": "fixed_reference_frequencies", "dropout": 0.2, "missing": "condition_on_call_available", "type": "generative_str"}`

## Ranking

| Rank | Model | logL | ΔlogL | Weight |
|---:|---|---:|---:|---:|
| 1 | hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity | -342.077934 | 0.000000 | 99.914% |
| 1 | belmonte_v1_scenario_a_all_loci / belmonte_v1_scenario_a_all_loci / amarna_no_identity | -342.077934 | 0.000000 | 99.914% |
| 2 | phizackerley_v2_scenario_a_all_loci / phizackerley_v2_scenario_a_all_loci / amarna_no_identity | -349.546457 | -7.468524 | 0.057% |
| 3 | variant_e1_scenario_a_all_loci / variant_e1_scenario_a_all_loci / amarna_no_identity | -350.496106 | -8.418172 | 0.022% |
| 4 | dodson_v2_scenario_a_all_loci / dodson_v2_scenario_a_all_loci / amarna_no_identity | -351.822134 | -9.744200 | 0.006% |
| 5 | tawfik_v1_scenario_a_all_loci / tawfik_v1_scenario_a_all_loci / amarna_no_identity | -353.704919 | -11.626986 | 0.001% |

## Support Ratios (LR for sensitivity settings; BF for population-prior comparisons)

| Comparison | Ratio |
|---|---:|
| hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity vs. belmonte_v1_scenario_a_all_loci / belmonte_v1_scenario_a_all_loci / amarna_no_identity | 1 |
| hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity vs. phizackerley_v2_scenario_a_all_loci / phizackerley_v2_scenario_a_all_loci / amarna_no_identity | 1752.02 |
| hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity vs. variant_e1_scenario_a_all_loci / variant_e1_scenario_a_all_loci / amarna_no_identity | 4528.62 |
| hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity vs. dodson_v2_scenario_a_all_loci / dodson_v2_scenario_a_all_loci / amarna_no_identity | 17055 |
| hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity vs. tawfik_v1_scenario_a_all_loci / tawfik_v1_scenario_a_all_loci / amarna_no_identity | 112082 |

## Locus Details

### Rank 1: hawass_v2_scenario_a_all_loci / hawass_v2_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -33.744855 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -63.745940 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.210276 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.542577 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.049760 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -56.741505 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.126568 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -28.916452 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 1: belmonte_v1_scenario_a_all_loci / belmonte_v1_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -33.744855 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -63.745940 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.210276 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.542577 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.049760 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -56.741505 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.126568 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -28.916452 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 2: phizackerley_v2_scenario_a_all_loci / phizackerley_v2_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -34.725462 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.795341 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.610109 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.111218 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.604459 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -58.820921 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.405058 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.473890 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 3: variant_e1_scenario_a_all_loci / variant_e1_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -34.725462 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 592704 |
| Locus2 | -65.642788 | blue, cyan, green, other, red | 15 | max factor rows 47250 |
| Locus3 | -37.610109 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 592704 |
| Locus4 | -48.896272 | blue, cyan, green, other, red, yellow | 21 | max factor rows 185220 |
| Locus5 | -36.320799 | blue, cyan, green, other, red, yellow | 21 | max factor rows 185220 |
| Locus6 | -58.820921 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 592704 |
| Locus7 | -39.005865 | blue, green, other, red, yellow | 15 | max factor rows 47250 |
| Locus8 | -29.473890 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 592704 |

### Rank 4: dodson_v2_scenario_a_all_loci / dodson_v2_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.418519 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -65.825119 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -38.303247 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -48.896276 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -35.552869 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -58.715573 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.405604 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -29.704928 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |

### Rank 5: tawfik_v1_scenario_a_all_loci / tawfik_v1_scenario_a_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.418525 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.795324 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -38.303251 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.804356 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -36.297417 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -59.514062 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.405063 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -30.166921 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
