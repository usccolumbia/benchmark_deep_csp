# Deep generative CSP models are implicit template methods: benchmark data, code and results

Data, code and per-target results for the paper **"Deep Generative Crystal Structure Prediction: A Benchmark Study and a Controlled Test of Prototype Dependence"**
(Lai Wei, Rongzhi Dong, Ying Feng, Madeline Miklos, Jianjun Hu; University of South Carolina).

The paper compares the template-based method **TCSPv3** with twelve deep generative CSP models on **T180** (180 Materials Project targets screened against the MP-20 training and
validation data) and **T46**, decomposes the generative models' successes against the template baseline, and removes whole stoichiometric prototype families from the
training data of CrystalFlow and EquiCSP to test prototype dependence causally.

## What is in this repository

| path | contents |
|---|---|
| `data/testsets/v3.1/{T180,T46}/` | reference CIFs, `metadata.csv`, the formula queries given to every method (`queries.csv/json`), `records.json`, `template_exclusions.json`, `manifest.json`; `CHANGES.md`, `CIF_SHA256SUMS.json`, `release_manifest.json` |
| `data/audit/` | the MP-20 leakage audit (report, per-target overlaps, flagged records) and the replacement records of the audit |
| `results/` | per-target scores for the 13 methods on T180 and T46 (`t180_per_case_all_methods.csv`, `t46_per_case_all_methods.csv`), bootstrap intervals, strict-primitive sensitivity, the decomposition against TCSPv3 and the prototype-availability counts, the TCSPv3 result with ORB-v3 relaxation, and the case studies (`case_studies/`, with the CIFs behind Figures 7 and 8) |
| `ablation/` | the removed-structure IDs of each family-removal split (`removed_ids/`), the audit record with split hashes, per-target and summary results of all ten arms (eight family arms for CrystalFlow and EquiCSP, plus the EquiCSP AB3 and perovskite-framework arms), and the parent-template analysis (`ablation_degeneracy_relaxed.csv`, `ablation_parent_detail_all.csv`, `ba6al10_framework_matches.csv`) |
| `scripts/` | the scripts that produced the results, tables and figures. They contain absolute paths of the authors' machines; set the path constants at the top of each file before running |
| `checkpoints/README.md` | the sources of the published checkpoints; the retrained ablation checkpoints are not released |
| GitHub Release assets | `T180_top20_predictions.zip`, `T46_top20_predictions.zip`, `ablation_top20_predictions.zip` (top-20 CIFs per method and target, with `MANIFEST.csv`) |

Prediction files are named `rankNN_<formula>_<candidate id>.cif` inside `<method>/<material id>/`; use `MANIFEST.csv` for paths. Rank means the first exported candidate:
ORB-v3 relaxation and energy-per-atom ranking for the generative models, TCSPv3's native chemical score (no relaxation) for TCSPv3.

## Protocol

Each method is given only the primitive-cell formula with integer counts of the target and returns up to 20 candidates. A candidate matches when pymatgen `StructureMatcher`
(`ltol` 0.2, `stol` 0.3, `angle_tol` 5, primitive cell, scaled, no supercell) finds it equivalent to the reference. Targets for which a method returns no valid candidate count as misses.
Space groups use `symprec` 0.1; consensus means the same top-1 candidate satisfies both criteria.

## T180 results (structure match; counts out of 180)

| Method | Top-1 | Top-5 | Top-20 |
|---|---:|---:|---:|
| TCSPv3 | 122 (67.8%) | 154 (85.6%) | 163 (90.6%) |
| *TCSPv3 + ORB-v3 relax* | *130 (72.2%)* | *159 (88.3%)* | *163 (90.6%)* |
| CrystalFlow | 115 (63.9%) | 121 (67.2%) | 130 (72.2%) |
| OMatG | 107 (59.4%) | 116 (64.4%) | 124 (68.9%) |
| Uni-3DAR | 105 (58.3%) | 110 (61.1%) | 114 (63.3%) |
| EquiCSP | 102 (56.7%) | 111 (61.7%) | 119 (66.1%) |
| DiffCSP | 97 (53.9%) | 108 (60.0%) | 115 (63.9%) |
| MatterGen | 97 (53.9%) | 110 (61.1%) | 119 (66.1%) |
| cond-CDVAE | 72 (40.0%) | 82 (45.6%) | 86 (47.8%) |
| GemsDiff | 62 (34.4%) | 64 (35.6%) | 66 (36.7%) |
| TGDMat | 50 (27.8%) | 52 (28.9%) | 53 (29.4%) |
| CDVAE | 48 (26.7%) | 52 (28.9%) | 53 (29.4%) |
| SymmCD | 40 (22.2%) | 45 (25.0%) | 49 (27.2%) |
| CrystalGRW | 23 (12.8%) | 24 (13.3%) | 24 (13.3%) |

The italic row relaxes and ranks TCSPv3's candidates with ORB-v3, the same post-processing the generative models receive; the plain TCSPv3 row uses its native ranking.

## T46 results (counts out of 46)

| Method | Top-1 | Top-5 | Top-20 |
|---|---:|---:|---:|
| TCSPv3 | 34 (73.9%) | 39 (84.8%) | 41 (89.1%) |
| Uni-3DAR | 32 (69.6%) | 32 (69.6%) | 33 (71.7%) |
| OMatG | 27 (58.7%) | 28 (60.9%) | 30 (65.2%) |
| EquiCSP | 27 (58.7%) | 29 (63.0%) | 29 (63.0%) |
| CrystalFlow | 26 (56.5%) | 27 (58.7%) | 31 (67.4%) |
| DiffCSP | 26 (56.5%) | 28 (60.9%) | 31 (67.4%) |
| MatterGen | 25 (54.3%) | 25 (54.3%) | 27 (58.7%) |
| cond-CDVAE | 17 (37.0%) | 18 (39.1%) | 19 (41.3%) |
| TGDMat | 17 (37.0%) | 17 (37.0%) | 18 (39.1%) |
| GemsDiff | 15 (32.6%) | 15 (32.6%) | 16 (34.8%) |
| CDVAE | 12 (26.1%) | 12 (26.1%) | 12 (26.1%) |
| CrystalGRW | 8 (17.4%) | 8 (17.4%) | 8 (17.4%) |
| SymmCD | 7 (15.2%) | 7 (15.2%) | 8 (17.4%) |

## Prototype-family ablation (T180 targets of the removed family; top-20 recovery)

| Model | Removed family | Targets | Baseline fit | Ablated fit | Ablated consensus | Degraded |
|---|---|---:|---:|---:|---:|---:|
| crystalflow | 1:1 binaries (AB) | 19 | 16 | 15 | 14 | 1 |
| crystalflow | 1:2 binaries (AB2) | 23 | 20 | 12 | 12 | 8 |
| crystalflow | 1:1:2 ternaries (ABC2) | 24 | 21 | 13 | 13 | 9 |
| crystalflow | 1:1:2:6 quaternaries (ABC2D6) | 11 | 8 | 10 | 9 | 0 |
| equicsp | 1:1 binaries (AB) | 19 | 16 | 15 | 14 | 1 |
| equicsp | 1:2 binaries (AB2) | 23 | 17 | 14 | 14 | 3 |
| equicsp | 1:1:2 ternaries (ABC2) | 24 | 21 | 13 | 13 | 10 |
| equicsp | 1:1:2:6 quaternaries (ABC2D6) | 11 | 8 | 10 | 8 | 0 |
| equicsp | 1:3 binaries (AB3) | 7 | 3 | 2 | 2 | 1 |
| equicsp | perovskite frameworks and the ABC2D6 family (PEROV), tested on the 11 ABC2D6 targets | 11 | 8 | 7 | 3 | 2 |

AB3 and PEROV were run for EquiCSP only. PEROV removes the ABC2D6 family and every training/validation structure whose framework matches a perovskite-type ABC2D6 target (2,819 training and 897 validation structures); its consensus recovery falls from 8 to 3 of the 11 targets.

*Baseline fit*: the model trained on the full MP-20 split recovers the target within its top 20. *Ablated fit* / *consensus*: the model retrained without the family. *Degraded*: baseline fit and no ablated fit.

## Reproducing

```bash
pip install -r requirements.txt
```

1. Predictions for new methods go in `<method>/<material id>/` folders as in the archives; `scripts/new_t180_topk_sg_all.py` scores them (top-k match rank, space-group agreement).
2. `scripts/new_t180_results_tables.py` and `scripts/new_t180_figures.py` build the tables and figures; `scripts/bootstrap_ci_t180.py` the intervals.
3. The ablation splits are rebuilt from the original MP-20 `train.csv`/`val.csv` (CDVAE release, 27,136 / 9,047 structures) by removing the IDs in `ablation/removed_ids/`.
   The audit record `ablation/audit.json` holds the SHA-256 of every split file.

## Citation

See `CITATION.cff`. The TCSPv3 software (version 3.0.0 for the results of this paper) is at https://github.com/usccolumbia/TCSPV3.
