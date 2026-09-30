# primary Evidence-Set Audit

Overall status: **PASS**

This audit checks the primary generative STR run files. It verifies that each run uses the explicit `generative_str` observation model, that all model families within a scenario use the same scored STR profile set, that KV21B is not present in the scored primary evidence set, and that the maternal profile attached to the KV62 fetuses is the same KV21A profile across model families even when historical labels differ.

## scenario_a_primary_models_all_loci_generative_str

Status: **PASS**

Observation model:

```json
{
  "type": "generative_str",
  "dropout": 0.2,
  "background_error": 0.001,
  "background_source": "fixed_reference_frequencies",
  "missing": "condition_on_call_available",
  "background_frequency_file": "../../data/background_allele_frequencies.json"
}
```

| Model | Scored labels | Fetal mother | KV21B present? |
|---|---|---|---|
| Hawass_v2_scenario_a | Akhenaten, Amenhotep, Fetus1, Fetus2, KV21A, KV35YL, Thuya, Tiye, Tutankhamun, Yuya | KV21A (`002447d937fd`) | no |
| Belmonte_v1_scenario_a | Amenhotep, Ankhesenamun, Baketaten, Fetus1, Fetus2, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Dodson_v2_scenario_a | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Phizackerley_v2_scenario_a | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Tawfik_v1_scenario_a | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Variant_E1_scenario_a | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |

## scenario_b_primary_models_all_loci_generative_str

Status: **PASS**

Observation model:

```json
{
  "type": "generative_str",
  "dropout": 0.2,
  "background_error": 0.001,
  "background_source": "fixed_reference_frequencies",
  "missing": "condition_on_call_available",
  "background_frequency_file": "../../data/background_allele_frequencies.json"
}
```

| Model | Scored labels | Fetal mother | KV21B present? |
|---|---|---|---|
| Hawass_v2_scenario_b | Akhenaten, Amenhotep, Fetus1, Fetus2, KV21A, KV35YL, Thuya, Tiye, Tutankhamun, Yuya | KV21A (`002447d937fd`) | no |
| Belmonte_v1_scenario_b | Amenhotep, Ankhesenamun, Baketaten, Fetus1, Fetus2, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Dodson_v2_scenario_b | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Phizackerley_v2_scenario_b | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Tawfik_v1_scenario_b | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Variant_E1_scenario_b | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |

## scenario_c_primary_models_all_loci_generative_str

Status: **PASS**

Observation model:

```json
{
  "type": "generative_str",
  "dropout": 0.2,
  "background_error": 0.001,
  "background_source": "fixed_reference_frequencies",
  "missing": "condition_on_call_available",
  "background_frequency_file": "../../data/background_allele_frequencies.json"
}
```

| Model | Scored labels | Fetal mother | KV21B present? |
|---|---|---|---|
| Hawass_v2_scenario_c | Akhenaten, Amenhotep, Fetus1, Fetus2, KV21A, KV35YL, Thuya, Tiye, Tutankhamun, Yuya | KV21A (`002447d937fd`) | no |
| Belmonte_v1_scenario_c | Amenhotep, Ankhesenamun, Baketaten, Fetus1, Fetus2, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Dodson_v2_scenario_c | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Phizackerley_v2_scenario_c | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Tawfik_v1_scenario_c | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Variant_E1_scenario_c | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |

## scenario_d_primary_models_all_loci_generative_str

Status: **PASS**

Observation model:

```json
{
  "type": "generative_str",
  "dropout": 0.2,
  "background_error": 0.001,
  "background_source": "fixed_reference_frequencies",
  "missing": "condition_on_call_available",
  "background_frequency_file": "../../data/background_allele_frequencies.json"
}
```

| Model | Scored labels | Fetal mother | KV21B present? |
|---|---|---|---|
| Hawass_v2_scenario_d | Akhenaten, Amenhotep, Fetus1, Fetus2, KV21A, KV35YL, Thuya, Tiye, Tutankhamun, Yuya | KV21A (`002447d937fd`) | no |
| Belmonte_v1_scenario_d | Amenhotep, Ankhesenamun, Baketaten, Fetus1, Fetus2, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Dodson_v2_scenario_d | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Phizackerley_v2_scenario_d | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Tawfik_v1_scenario_d | Amenhotep, Ankhesenamun, Fetus1, Fetus2, Meritaten, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
| Variant_E1_scenario_d | Amenhotep, Ankhesenamun, Fetus1, Fetus2, KV35YL, Smenkhkare, Thuya, Tiye, Tutankhamun, Yuya | Ankhesenamun (`002447d937fd`) | no |
