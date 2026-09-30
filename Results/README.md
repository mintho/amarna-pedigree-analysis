# Current calculation results

These are the retained numerical outputs for the documented analysis. Reproduction writes
to the separate, ignored `reproduced/` directory, not to this directory.

| Question | Start here |
|---|---|
| Which scored structure ranks highest? | [Structural ranking CSV](Primary_All_Locus/summary/structural_ranking_summary.csv) |
| How large are the support ratios? | [Primary summary](Primary_All_Locus/summary/report.md) |
| Which markers contribute to the gap? | [Per-locus CSV](Primary_All_Locus/summary/per_locus_vs_nearest.csv) |
| Does the full ranking survive perturbation and locus omission? | [All-locus sensitivity](Primary_All_Locus_Sensitivity/report.md) |
| Are local relationship rankings stable? | [Focused stability](Focused_Robustness/stability_report.md) |
| What changes when L1 or L2 is omitted? | [Focused jackknife](Focused_Robustness/jackknife_report.md) |
| How were the reference frequencies checked? | [Frequency audit](Frequency_prior_audit/audit.md) |
| Are the independent exact implementations consistent? | [Engine agreement](Validation/engine_agreement.json) |
| Where are the machine-readable numerical findings? | [Numerical facts](Manuscript_Tables/numerical_facts.json) |

All raw primary weights are structural weights. Historical aliases repeat the
same weight; do not add those aliases. A focused relationship weight compares
two local candidates and is not a five-structure historical model probability.
