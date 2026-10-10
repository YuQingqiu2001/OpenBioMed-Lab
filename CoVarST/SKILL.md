---
name: covarst
description: Run CoVarST fixed cancer references, H&E-to-spot inference and same-section reference-free spot-to-cell reconstruction.
---

# CoVarST operating contract

Read README.md and the relevant docs/ method/input schema before execution. Inspect `models/manifest.json`, run `verify-assets`, and select the exact cancer key. Never mix a checkpoint with another Program/gene order or reconstruct the fixed reference during inference.

Use an existing compatible environment. Installation or download requires the human user's authorization. This skill itself does not authorize messaging, publication, external uploads or dependency installation. Infer-spots accepts image features and coordinates only. Third-party weights are user-supplied and follow their official licenses.

For coarse ST reconstruction, use only current-section coarse RNA counts and real registered image nuclei. Never load external scRNA-seq, existing inferred cell expression or fine-resolution evaluation RNA into fitting/model selection. Preserve physical geometry, broad image class order, unassigned coarse mass and all-gene probability normalization. Continuous Program ranks are not cell types.

Use fresh output paths, record input hashes, and inspect complete markers/logs before claiming completion. Report student ensemble fidelity separately from the original LOO teacher scores. Conservation is an accounting result; biological cell-expression accuracy requires independent validation.

Commands and schemas: docs/QUICKSTART.md, docs/HE_TO_SPOTS.md, docs/REFERENCE_PROGRAMS.md and docs/SPOT_TO_CELL.md.
