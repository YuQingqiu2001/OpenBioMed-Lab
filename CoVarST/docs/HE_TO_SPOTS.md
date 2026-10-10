# H&E-only inference

The final cancer mapper retains the original architecture: 2560-D image input, hidden width 256, local and contextual morphology graphs (6 and 18 neighbors), GATv2 processing, RNA/composition heads and a learned support gate. There is one consolidated mapper per cancer. The student absorbs the teacher ensemble's accepted bounded INR effects into its targets and has no separate deployable INR network.

`extract-wsi` uses the original tissue-mask and physical hexagonal-grid mathematics: 112 μm patch field, 100 μm nearest-neighbor pitch and at least 0.20 tissue fraction with a tissue-positive center. It reads RGB patches, resizes to 224×224, applies the supplied ParamNet model and computes Virchow2 CLS + mean patch tokens (excluding its four register tokens), giving 2560 dimensions.

Feature NPZ schema:

| Key | Shape/type | Meaning |
|---|---|---|
| features | N × 2560, float16/float32 | ParamNet/Virchow2 features in the training representation |
| coords | N × 2, float32 | Consistent slide Cartesian coordinates |
| barcode | N unique strings | Exact feature/spot row order |
| coords_um | optional N × 2 | Explicit physical positions for WSI exports |

`infer-spots` loads no sample RNA, local theta or measured expression. It verifies checkpoint/reference hashes, loads the model state strictly, rebuilds the frozen graph, predicts Program probabilities and streams all 10,000 genes through `theta_rna @ W`. The neutral reference is the default; an explicit registered source×technology key can select a frozen batch W. Unknown keys fail rather than silently selecting an unrelated adapter.

Output H5 includes aligned coordinates/barcodes, theta RNA, theta composition, support probabilities, fixed gene IDs, expression probabilities and log1p(10,000×probability), with `transcriptomic_inputs_loaded=False` and complete/hash metadata. Probabilities are estimates of expression composition and do not recover a measured library size.

WSI extraction requires third-party weights obtained by the user from official sources. The public package does not include ParamNet, Virchow2 or CellViT++ weights and does not download them automatically. CellViT++ is needed for the separate real nucleus census in the spot-to-cell stage, not for decoding the cancer mapper's spot predictions.
