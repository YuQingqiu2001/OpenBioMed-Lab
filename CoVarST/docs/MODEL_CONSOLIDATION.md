# One unchanged-capacity checkpoint per cancer

The original final study contains 166 leave-one-patient-out mapper checkpoints. The release consolidates these into nine deployable cancer-specific students. It does not reduce hidden width, quantize weights, average raw network weights or bundle the old models behind a hidden ensemble.

For each same-cancer slide, every original teacher runs H&E-only inference in the same fixed Program coordinate system. Its accepted bounded INR correction is applied when present. RNA and composition softmax probabilities are averaged equally across all cancer folds. Averaging probabilities is valid here because all teachers share the same frozen W, Program order and 10,000-gene order.

One student with the original mapper capacity is initialized from a declared teacher checkpoint and trained on these frozen probability targets. The objective matches RNA/composition distributions, per-Program centered variation and local spatial differences. Every eligible slide and spot enters consolidation; features and geometry are exact barcode-aligned study inputs. Measured RNA is not newly read by the distillation pipeline, although teachers and fixed references were originally trained using RNA.

The initial teacher-target cache covers 302 eligible slides, 166 patients and 498,060 spots. One of the historical 303 slides was excluded by the original feature-coverage boundary. No scientific training-spot sampling is introduced by packaging. The teacher checkpoints and feature matrices remain local and are not included in the public package.

## Fidelity gates

The release requires patient/source-balanced mean RNA JSD ≤0.05, composition JSD ≤0.05, decoded-expression gene PCC ≥0.90, decoded spot cosine ≥0.98 and spatial theta-gradient cosine ≥0.85. All 10,000 genes enter decoded fidelity evaluation. Metrics and their balancing unit are recorded per slide and per cancer. Candidate selection enforces all gates, including the spatial gate, before preferring the lowest aggregate loss.

These are engineering consolidation gates measured against the ensemble on deployment consolidation inputs. Many teachers saw a given patient's RNA during their original training. The gates do not constitute a new held-out biological evaluation, do not prove noninferiority to original teachers on new patients and do not authorize reusing original LOO scores as student results.

Original metrics remain in `provenance/original_teacher_LOO_metrics.csv`; student results remain in `validation/` and per-cancer `model_card.json`. A model's Program/gene/reference hashes, architecture parameters, teacher fold count, training/selection contract and fidelity results accompany its checkpoint. A failed candidate is kept outside `models/`.

Independent patient-cohort validation of each consolidated student is still needed before assigning it biological or clinical performance claims.

## Reproduce consolidation locally

The optional portable recipe expects the original frozen study layout with its registry, all original teacher checkpoints, completed patient-out inference/spot-index files and registered H&E feature files. These large source inputs are not redistributed. Outputs go to a separate fresh workspace.

```bash
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase cache
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase targets
python scripts/consolidate_checkpoints.py --study-root ORIGINAL_STUDY --work-root WORK/consolidation --phase student
```

Run `targets` to completion before `student` when reproducibility and bounded resource sharing matter. The optional recipe preserves numerical losses and original architecture, relocates author paths, and retains a candidate satisfying every fidelity gate. Randomness/device arithmetic can prevent byte-identical checkpoints. Release weights remain identified by their actual recorded hashes and validation results.
