# Local release status

Nine cancer-specific consolidated checkpoints are complete and saved locally. Each cancer has exactly one `model_checkpoint.pt`, its matching Program group index and a model card. All nine satisfy every declared consolidation fidelity gate. Total mapper weight size is 189,067,520 bytes (189.1 MB); each mapper retains its original network capacity.

No remote push or GitHub release has been performed. The original 166 teacher checkpoints and all third-party weights remain outside this public subproject. Nine fixed references, the three workflow entry points, portable consolidation recipes and Chinese/English documentation are included.

## Consolidation fidelity

These patient/source-balanced statistics compare students with their same-cancer teacher ensemble on deployment consolidation inputs. They are not independent biological or clinical performance measurements. Original LOO teacher scores remain separately recorded in provenance/.

| Cancer | K | Teachers | Epoch | RNA JSD | Gene PCC | Spatial gradient cosine |
|---|---:|---:|---:|---:|---:|---:|
| bladder_cancer | 48 | 5 | 100 | 0.0008 | 0.9634 | 0.8512 |
| breast_cancer | 44 | 59 | 5 | 0.0052 | 0.9811 | 0.9318 |
| colorectal_cancer | 44 | 46 | 5 | 0.0121 | 0.9755 | 0.9003 |
| cutaneous_squamous_cell_carcinoma | 36 | 4 | 5 | 0.0167 | 0.9709 | 0.9187 |
| ependymoma | 32 | 11 | 5 | 0.0105 | 0.9660 | 0.8579 |
| kidney_clear_cell_carcinoma | 48 | 24 | 5 | 0.0160 | 0.9680 | 0.8665 |
| lung_adenocarcinoma | 48 | 5 | 15 | 0.0093 | 0.9589 | 0.8646 |
| pancreatic_adenocarcinoma | 48 | 4 | 60 | 0.0024 | 0.9691 | 0.8550 |
| prostate_adenocarcinoma | 48 | 8 | 25 | 0.0072 | 0.9737 | 0.8620 |

Full values: [student_fidelity_summary.csv](../validation/student_fidelity_summary.csv). Per-slide results and input/teacher/reference hashes accompany the package. Consolidation used all 302 eligible slides, 166 patients and 498,060 spots.

## Verification scope

Software checks cover engine imports, generic-coarse input packing, rejection of external RNA references, all-gene probability normalization, conserved-mass accounting and the strict nullspace projector against an independent SVD oracle. Real-data checks run each final checkpoint on one complete original feature slide, strictly load model state and verify all 10,000 decoded genes and exact gene/reference order. See [release_verification.json](../validation/release_verification.json) and [software_checks.json](../validation/software_checks.json).

Raw WSI extraction reuses the frozen tissue-grid, ParamNet and Virchow2 functions; the packaging validation did not rerun full WSI backbone extraction or retrain reference construction from raw BayesTME counts. The arbitrary-coarse adapter has software/input-contract validation and does not inherit independent biological acceptance. Independent patient-cohort and cell-expression validation remain outstanding.
