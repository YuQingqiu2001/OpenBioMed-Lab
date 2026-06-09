#!/usr/bin/env Rscript
# ============================================================
# export_bioc_genes.R — 从 Bioconductor 导出人+鼠基因符号
# 
# 用法:
#   Rscript export_bioc_genes.R [workspace]
#   
# 输出:
#   <workspace>/data/bioc_genes.json     (合并基因符号)
#   <workspace>/data/kegg_pathways.json  (KEGG 通路 + GO Biological Process)
#
# 依赖: org.Hs.eg.db, org.Mm.eg.db, jsonlite
# 可选: KEGGREST, GO.db, AnnotationDbi
# 环境: WSL2 conda sc_spatial_env
# ============================================================

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(jsonlite)
})

# --- 输出路径 ---
args <- commandArgs(trailingOnly = TRUE)
workspace <- if (length(args) >= 1) {
  args[[1]]
} else {
  Sys.getenv("LITERATURE_KG_ROOT", unset = file.path(getwd(), "literature-workspace"))
}
OUT_DIR <- file.path(workspace, "data")
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# --- 1. 人类基因符号 ---
cat("Extracting human gene symbols...\n")
hs_syms <- keys(org.Hs.eg.db, keytype = "SYMBOL")
# ALLCAPS 模式: A1BG, TNF, EGFR, BCL2, CD4
hs_pattern <- grep("^[A-Z][A-Z0-9]{1,7}$", hs_syms, value = TRUE, perl = TRUE)
cat(sprintf("  Human total: %d, ALLCAPS: %d\n", length(hs_syms), length(hs_pattern)))

# --- 2. 小鼠基因符号 ---
cat("Extracting mouse gene symbols...\n")
mm_syms <- keys(org.Mm.eg.db, keytype = "SYMBOL")
# TitleCase 模式: Tnf, Il6, Cd4, Trp53
mm_titlecase <- grep("^[A-Z][a-z0-9]{1,7}$", mm_syms, value = TRUE, perl = TRUE)
# ALLCAPS 模式 (少数): C2, C3, F5, F7
mm_allcaps   <- grep("^[A-Z][A-Z0-9]{1,7}$", mm_syms, value = TRUE, perl = TRUE)
mm_all <- unique(c(mm_titlecase, mm_allcaps))
cat(sprintf("  Mouse total: %d, TitleCase: %d, ALLCAPS: %d, Combined: %d\n",
            length(mm_syms), length(mm_titlecase), length(mm_allcaps), length(mm_all)))

# --- 3. 合并去重 ---
merged <- unique(c(hs_pattern, mm_all))
cat(sprintf("  Merged unique: %d\n", length(merged)))

# --- 4. 写入 JSON ---
write_json(merged, file.path(OUT_DIR, "bioc_genes.json"),
           auto_unbox = FALSE, pretty = FALSE)
cat(sprintf("Wrote %s (%d genes)\n",
            file.path(OUT_DIR, "bioc_genes.json"), length(merged)))

# --- 5. KEGG + GO Biological Process 术语 ---
cat("Extracting pathway and process terms...\n")
pathway_terms <- character()

if(requireNamespace("KEGGREST", quietly = TRUE)) {
  library(KEGGREST)
  hsa_paths <- keggList("pathway", "hsa")
  path_names <- gsub(" - Homo sapiens.*", "", hsa_paths)

  # 过滤疾病通路
  disease_terms <- c('cancer','carcinoma','tumor','tumour','melanoma','tuberculosis',
                     'leukemia','lymphoma','sarcoma','glioma','diabetes','alzheimer',
                     'parkinson','obesity','asthma','infection','disease','virus',
                     'hepatitis','influenza','measles','malaria','cholera','leishmania',
                     'toxoplasmosis','amebiasis','shigellosis','salmonella',
                     'pathogenic','viral')
  clean_kegg <- path_names[!grepl(paste(disease_terms, collapse = "|"),
                                   tolower(path_names))]
  pathway_terms <- c(pathway_terms, unname(clean_kegg))
  cat(sprintf("  KEGG total: %d, disease-filtered: %d\n",
              length(path_names), length(clean_kegg)))
} else {
  cat("  KEGGREST not available, skipping pathway export\n")
}

if(requireNamespace("GO.db", quietly = TRUE) &&
   requireNamespace("AnnotationDbi", quietly = TRUE)) {
  go_rows <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
    columns = c("TERM", "ONTOLOGY"),
    keytype = "GOID"
  )
  go_bp <- unique(go_rows$TERM[go_rows$ONTOLOGY == "BP" & !is.na(go_rows$TERM)])
  pathway_terms <- c(pathway_terms, go_bp)
  cat(sprintf("  GO Biological Process terms: %d\n", length(go_bp)))
} else {
  cat("  GO.db/AnnotationDbi not available, skipping GO term export\n")
}

pathway_terms <- unique(pathway_terms[nzchar(pathway_terms)])
write_json(pathway_terms, file.path(OUT_DIR, "kegg_pathways.json"),
           auto_unbox = FALSE, pretty = FALSE)
cat(sprintf("Wrote %s (%d pathway/process terms)\n",
            file.path(OUT_DIR, "kegg_pathways.json"), length(pathway_terms)))

cat("\nDone! Run 'python scripts/gen_edges.py' to regenerate edges.\n")
