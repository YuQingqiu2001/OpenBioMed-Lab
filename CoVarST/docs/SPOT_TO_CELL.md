# Reference-free spot-to-cell reconstruction

This stage uses the current section's coarse measured ST as its only RNA source. It also needs real image nuclei, their physical capture geometry, image features, broad PanNuke classes and positive morphology-derived capacities. No external scRNA-seq matrix, external cell-expression profile, previous cell expression or fine-resolution expression enters fitting. Broad lineage marker names provide weak anchors; their values are estimated from the same coarse ST.

The engine fits coarse-only class RNA profiles, shrinkage and 0–8 continuous axes per supported image class. Those axes are not cell types and are not the upstream cancer's fixed 32–48 Programs. Spatial blocked CV, buffers, image-feature cross fitting and permutation controls regulate the fit. The original cell refinement uses graph contrasts and a verified component measurement-nullspace projector. The absolute 1e-9 acceptance tolerance and all original measurement-row checks are retained.

Weakly supported classes remain rank zero under the original support rule. Global RNA profiles are estimated from this section; subsequent spatial CV is conditional/transductive and does not validate independent single-cell expression truth.

## Generic coarse input H5

The input adapter requires schema `covarst_coarse_counts_and_image_nuclei_v1`, `complete=1`, `rna_source="measured_coarse_st"`, `external_single_cell_reference_used=False`, and `real_image_nuclei=True`. It validates positive measured integer counts, unique IDs, physical coordinates and same-component geometry links.

| Dataset | Required shape / meaning |
|---|---|
| gene_name | G unique gene symbols in matrix order |
| matrix/{shape,data,indices,indptr} | CSC G × S measured coarse counts |
| geometry/{shape,data,indices,indptr} | CSR S × N real nucleus/Voronoi capture overlap weights |
| spots/barcode | S unique spot IDs |
| spots/coords_um | S × 2 physical Cartesian positions |
| spots/component | S nonnegative measured-domain component IDs |
| cells/cell_id | N unique integer nucleus IDs |
| cells/class_id | N integers 0–4: Neoplastic, Inflammatory, Connective, Dead, Epithelial |
| cells/coords_um | N × 2 physical Cartesian positions |
| cells/spatial_px | N × 2 registered image coordinates |
| cells/features | N × D standardized real nucleus image features |
| cells/rna_capacity | N positive morphology-derived capacities, never external RNA |
| cells/owner_bin | N spot indices for within-component ownership/interpolation |
| cells/domain_component | N component IDs matching owner spots |

The matrix columns, geometry rows and spot metadata must share order; geometry columns and cell metadata must share order. Preserve real overlap geometry, including orphan positive spots. A nearest-spot assignment does not substitute for capture-overlap geometry.

The generic adapter only packs inputs and independently estimates coarse profiles. It delegates numerical fitting to the frozen BayesGraph functions, with the same verified projector. Native16 and virtual55 preparation adapters preserve the historical physical affine, finite Voronoi domain and fractional virtual capture integration.

## Output semantics

Factor H5 contains class gene probabilities, all-gene loadings and cell state. The primary expression distribution is

`P_post = softmax_all_genes(log(class_gene_probability) + cell_state @ class_program_loadings)`.

`materialize-cells` computes this in blocks with a full-gene normalizer; a selected gene panel is never renormalized as if it were the full transcriptome. `log1p_cp10k` is an analysis representation of inferred probabilities.

`allocate-cell-mass` applies the frozen responsibility/IPF allocator to each original spot-gene observation. It writes fractional captured counts, spot-cell-gene links and unassigned gene mass for spots without nuclei. Its conservation checks establish accounting consistency. This output retains the coarse capture imprint and does not prove biological single-cell resolution; use `P_post` as the primary expression structure representation.

## Validation limits

Historical controlled HD calculations and image diagnostics are separate from biological acceptance. This release does not claim independently validated cell identity accuracy across cancer types. Fine 2 μm expression, when available, belongs to separate evaluation and must remain inaccessible to preparation/model selection. The generic arbitrary-coarse packaging adapter has software and input-contract checks; it is a packaging extension and does not inherit the original HD cohort's biological validation.
