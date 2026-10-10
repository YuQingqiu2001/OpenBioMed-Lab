# Third-party software and weights

Obtain third-party source and checkpoints from their official projects. No copies of their weights are included here and no automatic downloads occur.

| Component | Role | Official source |
|---|---|---|
| BayesTME 1.0.0 | Independent raw-count local Program decomposition | https://github.com/tansey-lab/bayestme |
| ParamNet | H&E stain normalization | https://github.com/khtao/ParamNet |
| Virchow2 | CLS + mean patch-token image features | https://huggingface.co/paige-ai/Virchow2 |
| CellViT++ | Real image nucleus segmentation and features | https://github.com/TIO-IKIM/CellViT-plus-plus |
| torch-geometric | GATv2 graph operations | https://github.com/pyg-team/pytorch_geometric |

Virchow2 access and use follow its model card/license. This package's license does not grant rights to any third-party checkpoint. ParamNet input root must contain the original `source.model.ParamNet` implementation and `checkpoints/ParamNet-Uni.pt`. Virchow2 root must contain its official compatible state. Follow CellViT++ documentation to produce a genuinely registered nucleus census; its class conventions must be mapped explicitly to the input schema.

The copied engine files are this study's original mathematical implementation/adapters. `provenance/source_manifest.json` records exact source hashes and where release packaging retains exact AST subsets or relocates runtime defaults. Historical source names are preserved internally for model compatibility.
