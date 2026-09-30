# Reference-frequency sources and retrieval

The calculations use modern reference-population STR allele frequencies from
the Huckenbeck-Scheil database at the University of Duesseldorf. These are
proxy population assumptions, not frequencies estimated from the mummies.
The prepared inputs and extraction results are included; complete third-party
webpages are not distributed.

## Source pages and population selection

The recorded URLs request the Wayback snapshot `20120920010159`
(20 September 2012, 01:01:59 UTC). An archive service may redirect to another
available capture. The JSON provenance records the requested source URLs.

| Label | Marker and source page | Selected population | Sample size |
|---|---|---|---:|
| L1 | [D13S317](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/d13s317.htm) | Egypt (pooled) | 260 |
| L2 | [D7S820](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/D7S820.html) | Egypt (pooled) | 260 |
| L3 | [D2S1338](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/D2S1338.html) | Israel (Jews) | 163 |
| L4 | [D21S11](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/d21s11.html) | Egypt (Cairo area) | 140 |
| L5 | [D16S539](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/D16S539.html) | Egypt (Central, El-Minia) | 120 |
| L6 | [D18S51](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/D18S51.html) | Egypt (Cairo area) | 140 |
| L7 | [CSF1PO](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/csf.html) | Egypt (pooled) | 219 |
| L8 | [FGA](https://web.archive.org/web/20120920010159/http://www.uni-duesseldorf.de/WWW/MedFak/Serology/DNA-Systeme/fga.htm) | Egypt (pooled) | 390 |

The L3 fallback is explicit because the selected database page did not supply
an Egyptian D2S1338 frequency table. Original allele labels, selected
populations, sample sizes and available reference identifiers are retained in
[the extracted frequency dataset](../data/allele_frequencies.duesseldorf_egypt.json).
The [audit](../../Results/Frequency_prior_audit/audit.md) and
[mapping report](../../Results/Frequency_prior_audit/color_mapping.md) explain
how these values relate to the operational locus-specific allele tokens.

## Re-extract the sources

From the repository root, run:

```bash
python3 tools/reproduce.py priors
```

The extractor downloads the eight pages if they are not already present in
`reproduced/source_cache/duesseldorf_archive/`. That directory is ignored by
Git. Subsequent invocations reuse the locally retrieved pages. Extracted JSON
and audit reports are written to `reproduced/Frequency_prior_audit/`; the
retained inputs in `Data/` are not overwritten.

This optional provenance audit requires network access on first retrieval and
depends on the archive's availability. It does not authenticate the DNA or
automatically accept new source data. Review extracted frequencies and
population selection before using any changed inputs. Archive-injected HTML
can vary between downloads, so compare parsed frequency values and provenance,
not only webpage bytes. A parsing failure or unavailable source must be
reported rather than silently replaced with a different population.

All numerical calculations can be reproduced offline from the included JSON:

```bash
python3 tools/reproduce.py check
python3 tools/reproduce.py summary
python3 tools/reproduce.py primary
python3 tools/reproduce.py focused
python3 tools/reproduce.py sensitivity
```

`python3 tools/reproduce.py all` also invokes the source-retrieval audit and
therefore needs either a local ignored cache or archive access.

## Source-copy fingerprints

These SHA-256 digests identify the local HTML copies used for the retained
extraction. The copies themselves are not included. A new Wayback response
can differ in archive-added markup without changing the scientific table.

| Original page | SHA-256 |
|---|---|
| `d13s317.htm` | `c20687e71bb4eb1487bc21d52ec6419dc989959edfe2d31fab1e15501520e562` |
| `D7S820.html` | `df55340763880227ae8432d6d106d20ce2e30adb52ca41109a2c46b47c47c6d0` |
| `D2S1338.html` | `ded5ed252e5244fb1c3c3437d0212013f3657db3229dfcad95222ecbcab6c05b` |
| `d21s11.html` | `73d8444bb5dbed7e9848cdeaed80e70b8ce8a079b51a0ac737691a3953f2888c` |
| `D16S539.html` | `d51f6338a8a3e3c03f565aadf229cdf55d2aed93c0306820a6af85462bfdd245` |
| `D18S51.html` | `b6bc0b9c123dc3eb0f83b4c3017cef4a1e9fb26804512ade1e34bf347f73f40f` |
| `csf.html` | `a8bc64b599e5a10bf647208fddd988f4fb655d1bfaa4da460bd1c9a2b03b813d` |
| `fga.htm` | `3a9b9a9800467c37c76d7e2a7b6d5a44d7fb3ae68299d35fe52374c1b36e5d0c` |

The source pages carry Huckenbeck-Scheil copyright notices. Retrieval links
and scientific attribution are not a grant to redistribute the pages or a
licence to their underlying database. See [licensing scope](../../LICENSING.md).
