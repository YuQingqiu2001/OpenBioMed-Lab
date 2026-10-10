# Quick start

Run all commands inside the `CoVarST` checkout with `PYTHONPATH=src`. Inputs and outputs are explicit. Outputs require fresh paths. Dependencies are never installed automatically.

## H&E → spot ST

```bash
python -m covarst extract-wsi --wsi DATA/slide.svs \
  --paramnet-root THIRD_PARTY/ParamNet --virchow2-root THIRD_PARTY/Virchow2 \
  --output WORK/slide.features.npz --device cuda --batch-size 16
python -m covarst infer-spots --cancer colorectal_cancer \
  --features WORK/slide.features.npz --output WORK/slide.spot_expression.h5 --device cuda
```

MPP is read from the slide metadata, with the historical objective-power fallback recorded. Provide `--mpp` only when a measured physical scale is known. Existing 2560-D feature NPZs can bypass extraction. Use `--batch SOURCE|TECHNOLOGY` only for an exact declared fixed adapter; otherwise use the neutral W.

## Generic coarse ST → cells

Prepare a HDF5 using the schema in [SPOT_TO_CELL.md](SPOT_TO_CELL.md). A capture geometry matrix must link the same measured spots to real image nuclei; class IDs must use the declared five-class PanNuke order.

```bash
python -m covarst prepare-coarse --input DATA/section.coarse_and_nuclei.h5 \
  --sample SECTION --output WORK/section.prepared.h5
python -m covarst fit-coarse -- --sample SECTION --prepared WORK/section.prepared.h5 \
  --output-dir WORK/section.factors --ranks 0,8 --device cuda --max-cuda-gib 6
python -m covarst materialize-cells \
  --factors WORK/section.factors/SECTION.V55.BayesGraph_rank8_all_gene.factorized.h5 \
  --output WORK/section.P_post.h5
python -m covarst allocate-cell-mass --prepared WORK/section.prepared.h5 \
  --factors WORK/section.factors/SECTION.V55.BayesGraph_rank8_all_gene.factorized.h5 \
  --output WORK/section.X_mass.h5
```

The historical fitter retains `V55` in the factor filename for compatibility. The input schema and factor summary explicitly record whether this was genuine coarse ST or a virtual55 degradation. The `0,8` sweep returns separate rank candidates; choose by the recorded coarse spatial CV, never by fine-expression test scores. Biological cell identity is not established by this CV.

## Controlled HD degradation

```bash
python -m covarst prepare-virtual55 --help
python -m covarst fit-virtual55 --help
python -m covarst prepare-native16 --help
python -m covarst fit-native16 --help
```

`prepare-virtual55` consumes raw native16 counts, positions, scalefactors and real nucleus feature H5. Its virtual55 fractional capture mass is fitted independently. It has no fine-expression input. `prepare-native16 --mode coarse` fits a native16 census/profile package; the strict native16 fitter uses the memory-bounded rank-revealing projector. Fine-resolution data are reserved for separately invoked evaluation.

## Build a new cancer reference

```bash
python -m covarst prepare-reference-counts --help
python -m covarst deconvolve-local --help
python -m covarst fit-fixed-reference --help
python -m covarst refine-reference --help
python -m covarst train-he --help
```

Run the BayesTME command with a separately prepared BayesTME 1.0.0 environment. Explicitly specify the frozen candidate set and selection settings in [REFERENCE_PROGRAMS.md](REFERENCE_PROGRAMS.md); historical engine defaults belong to earlier development runs. A new cancer needs its own new reference and mapper training and cannot be assigned another cancer's gene/program order.
