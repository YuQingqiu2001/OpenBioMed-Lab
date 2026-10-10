# CoVarST

CoVarST converts H&E image features to spot-level spatial RNA using a fixed, cancer-specific Program reference. Its separate spot-to-cell module reconstructs cell-resolution expression from a section's own coarse ST counts and real image nuclei, without an external scRNA-seq reference.

This local release contains nine consolidated cancer-specific checkpoints, all passing the declared engineering fidelity gates. Check [release status](docs/RELEASE_STATUS.md) and `models/manifest.json` for exact assets and validation scope. The original 166 leave-one-patient-out checkpoints stay outside this package. No remote push has been performed.

| Public component | Contents |
|---|---|
| Reference construction | Raw-count local BayesTME decomposition; local-to-consensus mapping; candidate selection; fixed-reference joint refinement; nine frozen references and gene orders |
| H&E → spot ST | One final mapper checkpoint per cancer; ParamNet → Virchow2 feature extraction; dual-head graph mapper; frozen-reference all-gene decoding |
| Coarse ST → cells | Same-section RNA profiles, 0–8 continuous within-class axes, image geometry/features, spatial cross fitting, verified measurement-nullspace projection, all-gene probabilities and a separate conserved-mass ledger |

Cancer keys: `bladder_cancer`, `breast_cancer`, `colorectal_cancer`, `cutaneous_squamous_cell_carcinoma`, `ependymoma`, `kidney_clear_cell_carcinoma`, `lung_adenocarcinoma`, `pancreatic_adenocarcinoma`, `prostate_adenocarcinoma`.

[中文指南](GUIDE_ZH.md) · [Quick start](docs/QUICKSTART.md) · [Reference construction](docs/REFERENCE_PROGRAMS.md) · [HE inference](docs/HE_TO_SPOTS.md) · [Spot-to-cell](docs/SPOT_TO_CELL.md) · [Model consolidation and validation](docs/MODEL_CONSOLIDATION.md) · [Third-party dependencies](docs/THIRD_PARTY.md) · [Agent instructions](SKILL.md)

## Run from the checkout

Use an existing compatible environment, or create one explicitly. No dependency or backbone is downloaded automatically. The locally checked environment is Python 3.9, PyTorch 2.7.1, torch-geometric 2.6.1, numpy 1.26.4, scipy 1.10.1 and timm 0.9.16. BayesTME 1.0.0 runs separately in a compatible environment.

```bash
cd CoVarST
export PYTHONPATH="$PWD/src"
python -m covarst --help
python -m covarst list-models
python -m covarst verify-assets
python -m covarst infer-spots --cancer colorectal_cancer \
  --features DATA/slide.features.npz --output WORK/slide.spot_expression.h5 --device cuda
```

Installing this package with `pip install -e .` is optional. Keep the `models/` and `references/` directories together; use `--assets /path/to/CoVarST` when running elsewhere.

## Interpret the outputs

Spot expression is a decoded probability or log1p CP10K estimate on a fixed 10,000-gene panel. It is not measured UMI counts. Cell `P_post` is an inferred, all-gene-normalized expression distribution; `X_mass` is fractional captured-count allocation with an explicit unassigned-mass ledger. Count conservation does not establish true cell expression.

The original patient-out mapper scores in `provenance/` describe the original teachers under a fixed, cohort-derived reference. Consolidated student fidelity measures imitation of a same-cancer teacher ensemble and does not inherit those biological performance scores. Independent biological/clinical validation of the consolidated students remains outstanding.

No patient RNA matrices, H&E images, original teacher checkpoint collection or third-party weights are bundled. Frozen references and identifiers follow the historical study provenance. Historical `FSegGATv2` class names remain inside the compatibility engine to preserve checkpoint loading; the public method name is **CoVarST**.

## Layout

```text
configs/       runnable recipe configurations and broad lineage anchors
models/        one checkpoint + Program group index + model card per cancer
references/    frozen W, gene order, adapters, local-to-consensus maps
src/covarst/   public CLI and frozen mathematical engine
provenance/    source hashes, original teacher registry and original LOO scores
validation/    student fidelity and software verification reports
docs/          methods, schemas, execution and validation limits
```

Code is provided under [CC BY-NC-SA 4.0](LICENSE). Third-party components retain their own licenses and access conditions; obtain their files from the official links. Manuscript citation metadata will be added when available.
