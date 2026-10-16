# T180 stoichiometric-family ablation

Approved scope: CrystalFlow, Crystalite and EquiCSP, each with AB2, ABC2, AB and
ABC2D6 removal arms (12 new training runs). Reuse the existing Enigma pretrained
baseline checkpoints and their corrected T180 inference runs; do not retrain
baselines. No random-removal controls.
TCSPV3 stays unchanged as a template reference; whole-family deletion would
disable its same-stoichiometry retrieval by construction.

Families are sorted integer element counts divided by their greatest common
divisor: (1,2), (1,1,2), (1,1), (1,1,2,6). Element identities and their ordering do
not affect membership. Both train and validation are filtered, with the original
CSV fields retained and CIF-derived compositions audited. Test data are not
filtered and are never used for training or checkpoint selection.

All models must use the fixed original CDVAE MP20 release membership:
27136 training, 9047 validation, 9046 test. The unchanged DiffCSP CSV files
retain those exact split memberships and are used for their compatible feature
columns. Do not concatenate/re-split, move boundary rows, or adopt an
algorithm's reconstructed split. Release Git-blob hashes and membership
comparisons are recorded in original_mp20_verification.json. The approved
filtered datasets for all three models are data_selected/splits/.

The frozen T180 queries come from T180_campaign_20261006, including the resolved
20-atom Sr4Tm2Nb2O12 input for mp-557387. Inference uses explicit primitive-cell
counts, not reduced formulas. Four targets exceed the MP20 20-atom training
range and must be separately labeled as size extrapolation, not silently dropped.

Each algorithm/arm has separate caches and output directories. Fresh training
never loads a published checkpoint. Published checkpoint hparams specify the
CrystalFlow and EquiCSP architecture/optimizer configuration only. Both use
64-structure microbatches and accumulation 4 for effective batch 256, FP32,
seed 42. Fixed optimizer budgets: CrystalFlow 318000 and EquiCSP 106000
(3000/1000 baseline epochs, respectively, at 106 optimizer updates per epoch).
Filtered arms run to the same number of optimizer steps. Native validation-loss
LR scheduling is retained; validation-selected checkpoints are primary and
final-step checkpoints are secondary. Validation contains no removed-family structures.

Selected by the user on 2026-10-06: AB2 (23 T180 targets), ABC2 (24), AB (19),
ABC2D6 (11). These are four independent experiments, not one combined removal.
The original AB3/ABC4 preparation is retained as historical data, but is not an
approved training arm. The approved filtered splits live in data_selected/.

Crystalite retains its baseline 1024-width, 14-layer, 16-head architecture,
batch 128, seed 123 and 1400000-update checkpoint-matched budget. Its baseline
metadata uses a 5000000-update learning-rate schedule; retain that schedule and
cap execution at 1400000 updates. V100 uses FP32 rather than baseline BF16.
The reconstructed 27138/9046 split and data_crystalite/ are superseded and
retained only for diagnostic provenance; they must not supply manuscript
results. Crystalite's corrected runs use fresh caches/output directories
under training_original_mp20/crystalite/. No capacity reduction is allowed.
Reused author checkpoints do not prove complete original
training provenance; see BASELINE_REUSE.md for comparability limitations.

First submit actual V100 training smoke jobs, retaining model capacity. These
are feasibility checks only and never benchmark results. Full jobs follow only
after each algorithm's gate passes. All jobs use gpu-v100-16gb, one GPU, 4 CPUs.

Primary manuscript deliverables are Tables 3 and 4, matching the supplied table
definitions rather than reproducing historical numerical values from another
cohort. Table 3 uses one common fixed original MP20 split for all models.
Table 4 has four family-removal rows per generative algorithm (12 comparisons).
The fixed family denominators are AB2=23, ABC2=24, AB=19, ABC2D6=11.

For Table 4, predict only the corresponding removed-family subset for each
ablated model: 77 targets per algorithm, 231 target/model predictions total.
Reuse the same target IDs from the existing Enigma T180 baseline predictions.
Full-T180 ablated inference is optional and is not required for these tables.

Evaluation: generate 20 candidates per relevant target using the same audited
model-specific step counts and postprocessing in baseline/removal arms. Preserve
Crystalite's recorded sampling seed. Historic CrystalFlow/EquiCSP inference RNG
states are not established: do not infer them from training seeds or claim
seed-paired draws, and do not rerun baselines without authorization. See
baseline_adapter_audit_20261006/INFERENCE_CONTRACT.md for the inspected settings
and this comparability limitation. Report raw
StructureMatcher top-1/5/20 recovery, then separately any common relaxation and
ranking results. Joint recovery requires the same candidate to pass structural
matching and space-group comparison. Degraded counts baseline-success AND
ablated-failure on the removed family's fixed T180 subset. Missing generation is
reported explicitly. Report evaluated-family coverage alongside family metrics.
Table 4 ablated fit is the number of family targets matched by at least one of
the 20 candidates. Ablated consensus additionally requires matching inferred
space group for that same candidate. Degraded is baseline StructureMatcher
success AND ablated StructureMatcher failure, divided by the entire family's
fixed target count (not by the baseline-success count). Supplementary data
should retain baseline fit, gained/lost target IDs, candidate provenance and
missing/failed cases so every table cell can be independently recomputed.

Without random-removal controls, decreases cannot uniquely separate family
removal from the effect of having fewer distinct training examples. One seed
does not establish robustness across training randomness.

## Addendum (2026-10-15/16): EquiCSP AB3 and perovskite-framework arms

Two further EquiCSP arms were trained on Enigma with the same recipe as the earlier arms (106000 optimizer steps, effective batch 256,
seed 42, validation-selected checkpoint) and from fresh initialisation: AB3 (1:3 binaries; 1620 training and 518 validation structures removed; 7 T180 targets)
and PEROV (the ABC2D6 family plus every training/validation structure whose species-erased framework matches a perovskite-type ABC2D6 target;
2819 and 897 structures removed; evaluated on the 11 ABC2D6 targets). Inference: 20 candidates per target, ORB-v3 relaxation and energy ranking;
baseline: the full-data EquiCSP of the T180 campaign. CrystalFlow was not run for these two arms. Each arm is a single training run (seed 42).
