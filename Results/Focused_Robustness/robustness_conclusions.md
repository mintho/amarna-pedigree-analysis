# Robustness Conclusions

Stability report: `Results/Focused_Robustness/stability_report.json`

## Classification Rules

| Label | Rule |
|---|---|
| `stabil` | rank_flip_count == 0 and baseline_best_fraction == 1.0 |
| `annahmeabhaengig` | rank flips occur, but baseline_best_fraction >= 0.5 |
| `instabil` | baseline_best_fraction < 0.5 |

## Conclusions

| Run | Classification | Baseline best in matrix | Rank flips | Main driver |
|---|---|---:|---:|---|
| `focused_yuya_thuya_tiye_trio_l1_l2_run.json` | `stabil` | 100.0% | 0/27 | `keine` |
| `focused_amenhotep_tiye_children_l1_l2_run.json` | `stabil` | 100.0% | 0/27 | `keine` |
| `focused_yuya_kv55_grandparent_l1_l2_run.json` | `stabil` | 100.0% | 0/27 | `keine` |

## Interpretation

### `focused_yuya_thuya_tiye_trio_l1_l2_run.json`

`focused_yuya_thuya_tiye_trio_l1_l2_relationships / amarna_no_identity` bleibt in allen 27 Matrixpunkten Rang 1. Das Ergebnis ist innerhalb des gerechneten Annahmeraums stabil.

### `focused_amenhotep_tiye_children_l1_l2_run.json`

`focused_amenhotep_tiye_children_l1_l2_relationships / amarna_no_identity` bleibt in allen 27 Matrixpunkten Rang 1. Das Ergebnis ist innerhalb des gerechneten Annahmeraums stabil.

### `focused_yuya_kv55_grandparent_l1_l2_run.json`

`focused_yuya_kv55_grandparent_l1_l2_relationships / amarna_no_identity` bleibt in allen 27 Matrixpunkten Rang 1. Das Ergebnis ist innerhalb des gerechneten Annahmeraums stabil.
