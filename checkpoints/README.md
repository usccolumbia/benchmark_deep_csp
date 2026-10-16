# Checkpoints

## Ablation (family-removal) checkpoints

The retrained ablation checkpoints (CrystalFlow and EquiCSP, one per removed family, plus the EquiCSP AB3 and perovskite-framework models) are **not released**. The removed-structure IDs, split hashes, training settings and every generated candidate are provided, so the arms can be retrained from the original MP-20 split (`ablation/PROTOCOL.md`, `ablation/audit.json`).

## Other checkpoints used in the benchmark

Retrained or adapted by us (hosted on Figshare) and published checkpoints:

| Method | Checkpoint | Location |
|---|---|---|
| CrystalFlow | CSP-mp-20 | https://github.com/ixsluo/CrystalFlow/releases |
| OMatG | MP-20-CSP (Linear-ODE) | https://huggingface.co/OMatG |
| Uni-3DAR | mp20_csp.pt (released MP-20 CSP model) | https://huggingface.co/dptech/Uni-3DAR/tree/main |
| EquiCSP | pretrained mp_csp model | https://drive.google.com/drive/folders/1Dw945rCLCT3jbQ8FZsAxCCVMjeo-PByc?usp=drive_link |
| MatterGen | mattergen_run1, epoch 899, trained on the MP-20 training split | https://doi.org/10.6084/m9.figshare.34336494 |
| DiffCSP | pretrained mp_csp model | https://drive.google.com/drive/folders/11WOc9lTZN4hkIY7SKLCIrbsTMGy9TsoW?usp=sharing |
| cond-CDVAE | condcdvae-4m | https://github.com/ixsluo/cond-cdvae/releases/tag/v1.0.0 |
| GemsDiff | without_cell_diffusion model trained in February 2024 | https://doi.org/10.6084/m9.figshare.34336932 |
| TGDMat | MP-20 CSP-task checkpoint (model_400) | https://doi.org/10.6084/m9.figshare.34336773 |
| SymmCD | symmcd_run1 | https://doi.org/10.6084/m9.figshare.34336500 |
| CrystalGRW | pretrained_alexmp20 (Alex-MP-20); unconditional model with atom types and count forced to the query | https://zenodo.org/records/14948252 |
| CDVAE | trained on the MP-20 training split | https://doi.org/10.6084/m9.figshare.34336497 |
