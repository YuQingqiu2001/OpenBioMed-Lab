# From independent local Programs to a cancer reference

The final artifact is `references/CANCER/consensus_reference_probability_10k.npy`, W with shape K × 10,000. Each row is a gene probability distribution. `gene_order_10k.npy`, batch adapters and Program order are immutable model assets. Both the original teachers and the consolidated cancer student decode theta through that same W.

1. Pool measured raw-count slides by cancer. Build a cancer-specific 10,000-gene panel using marker coverage followed by equal-slide logCP10K-variance HVGs. Exclude mitochondrial, RPS and RPL genes as in the final configuration. Save exact original gene IDs/order and slide provenance; do not replace counts by normalized values for BayesTME.
2. Run unmodified upstream BayesTME 1.0.0 separately on each slide with local K in 2–12. Preserve the local posterior theta and local gene basis W. These local Programs are not shared identities across slides.
3. Fit soft local-to-consensus mappings and candidate shared dictionaries. The final cancer candidate set is 32, 36, 40, 44, 48, with preferred K=40. Select among candidates within 0.002 of the best validation gene PCC, preferring the size nearest 40. Preserve the candidate scores and selected mapping. The RCCD helper provides mapping machinery; the final locked-split direct-consensus entry point applies the recorded candidate-selection rule.
4. Jointly refine the selected shared W, free spot theta and bounded technology/source×technology adapters. The final objective combines logCP10K Huber loss, gene PCC, cosine and a weak theta anchor. This refinement is not accurately described as merely averaging local bases or as a new KL-NMF fit. Adapter variation is tied to platform/batch metadata, not patient identities.
5. Freeze K, W, Program order, gene order and adapter arrays with hashes. Train and evaluate image mappers against this reference without rebuilding it per fold.

| Cancer key | Fixed K |
|---|---:|
| bladder_cancer | 48 |
| breast_cancer | 44 |
| colorectal_cancer | 44 |
| cutaneous_squamous_cell_carcinoma | 36 |
| ependymoma | 32 |
| kidney_clear_cell_carcinoma | 48 |
| lung_adenocarcinoma | 48 |
| pancreatic_adenocarcinoma | 48 |
| prostate_adenocarcinoma | 48 |

The nine references total 396 Programs. The breast preparation followed its historical breast-specific input lineage; its frozen 44-Program assets are included. The generic preparation example covers the eight later cancer batches, and is an input template rather than a claim that all nine original data pipelines used identical paths.

## Executable recipes and records

`prepare-reference-counts`, `deconvolve-local`, `fit-fixed-reference`, `refine-reference` and `train-he` expose the frozen upstream CLIs. Their `--help` describes required measured-count queues, spot-index tables, split manifests and output directories. Use explicit candidates:

```bash
python -m covarst fit-fixed-reference --queue-dir WORK/local_bayes \
  --master-spot-index WORK/master_spot_index.csv --split-manifest WORK/reference_split.csv \
  --output-dir WORK/shared_reference --candidates 32 36 40 44 48 \
  --target-programs 40 --selection-tolerance 0.002
```

`compress-programs` exposes the earlier RCCD-only helper for reproducibility; its development defaults are not the final recipe. Use the locked-split direct-consensus recipe above for final candidate selection, then the recorded joint refinement. Exact per-cancer configuration and historical selection/refinement summaries accompany each frozen reference. `local_to_consensus_mapping.npy`, `local_program_assignment.csv`, `candidate_dictionary_selection.csv` and `consensus_contract.json` make the compression auditable. Patient theta targets are deliberately excluded from the public package.

## Evaluation scope

The historical fixed reference was constructed using the full available cancer RNA cohort before leave-one-patient-out image-mapper evaluation. The reference is therefore transductive. The original LOO results characterize image generalization conditional on that fixed reference and do not establish fully RNA-independent reference generalization. Consolidating teachers does not change this boundary.
