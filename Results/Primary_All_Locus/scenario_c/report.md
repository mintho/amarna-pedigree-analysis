# Pedigree Likelihood Report

Run configuration: `Data/models/runs/scenario_c_primary_models_all_loci_generative_str.json`
Engine: `factors`
Observation model: `{"background_error": 0.001, "background_frequency_file": "../../data/background_allele_frequencies.json", "background_source": "fixed_reference_frequencies", "dropout": 0.2, "missing": "condition_on_call_available", "type": "generative_str"}`

## Ranking

| Rank | Model | logL | ΔlogL | Weight |
|---:|---|---:|---:|---:|
| 1 | hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity | -348.363421 | 0.000000 | 99.979% |
| 1 | belmonte_v1_scenario_c_all_loci / belmonte_v1_scenario_c_all_loci / amarna_no_identity | -348.363421 | 0.000000 | 99.979% |
| 2 | phizackerley_v2_scenario_c_all_loci / phizackerley_v2_scenario_c_all_loci / amarna_no_identity | -356.899706 | -8.536285 | 0.020% |
| 3 | dodson_v2_scenario_c_all_loci / dodson_v2_scenario_c_all_loci / amarna_no_identity | -360.086805 | -11.723384 | 0.001% |
| 4 | variant_e1_scenario_c_all_loci / variant_e1_scenario_c_all_loci / amarna_no_identity | -360.995948 | -12.632527 | 0.000% |
| 5 | tawfik_v1_scenario_c_all_loci / tawfik_v1_scenario_c_all_loci / amarna_no_identity | -361.058234 | -12.694813 | 0.000% |

## Support Ratios (LR for sensitivity settings; BF for population-prior comparisons)

| Comparison | Ratio |
|---|---:|
| hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity vs. belmonte_v1_scenario_c_all_loci / belmonte_v1_scenario_c_all_loci / amarna_no_identity | 1 |
| hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity vs. phizackerley_v2_scenario_c_all_loci / phizackerley_v2_scenario_c_all_loci / amarna_no_identity | 5096.38 |
| hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity vs. dodson_v2_scenario_c_all_loci / dodson_v2_scenario_c_all_loci / amarna_no_identity | 123424 |
| hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity vs. variant_e1_scenario_c_all_loci / variant_e1_scenario_c_all_loci / amarna_no_identity | 306363 |
| hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity vs. tawfik_v1_scenario_c_all_loci / tawfik_v1_scenario_c_all_loci / amarna_no_identity | 326052 |

## Locus Details

### Rank 1: hawass_v2_scenario_c_all_loci / hawass_v2_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.996722 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.205479 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.402165 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.904980 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.748671 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -58.948697 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.007327 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.149380 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 1: belmonte_v1_scenario_c_all_loci / belmonte_v1_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -35.996722 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -64.205479 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -37.402165 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -47.904980 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -35.748671 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -58.948697 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.007327 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.149380 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 2: phizackerley_v2_scenario_c_all_loci / phizackerley_v2_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -36.689857 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -65.401231 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -38.369662 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.717606 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -36.353361 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -61.028173 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.513739 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.826077 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 3: dodson_v2_scenario_c_all_loci / dodson_v2_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -37.382896 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus2 | -66.766915 | blue, cyan, green, other, red | 15 | max factor rows 3375 |
| Locus3 | -39.062801 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus4 | -49.679640 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus5 | -36.441953 | blue, cyan, green, other, red, yellow | 21 | max factor rows 9261 |
| Locus6 | -60.922765 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |
| Locus7 | -39.513711 | blue, green, other, red, yellow | 15 | max factor rows 3375 |
| Locus8 | -30.316124 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 21952 |

### Rank 4: variant_e1_scenario_c_all_loci / variant_e1_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -36.689857 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -69.271894 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -38.369662 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -49.679661 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -37.153677 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -61.028173 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -38.976947 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -29.826077 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |

### Rank 5: tawfik_v1_scenario_c_all_loci / tawfik_v1_scenario_c_all_loci / amarna_no_identity

| Locus | logL | Alleles | Genotype states | Work summary |
|---|---:|---|---:|---:|
| Locus1 | -37.382932 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus2 | -65.401220 | blue, cyan, green, other, red | 15 | max factor rows 50625 |
| Locus3 | -39.062805 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus4 | -50.410748 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus5 | -37.046346 | blue, cyan, green, other, red, yellow | 21 | max factor rows 194481 |
| Locus6 | -61.721314 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
| Locus7 | -39.513757 | blue, green, other, red, yellow | 15 | max factor rows 50625 |
| Locus8 | -30.519110 | blue, cyan, green, other, pink, red, yellow | 28 | max factor rows 614656 |
