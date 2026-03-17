# ============================================================
# LRP1-associated PROGENy (pseudobulk) across:
#   1) HLCA fibro (fibro_subset)  : donor-level pseudobulk
#   2) Mouse BLM fibro (makoto.figure3): mouse-level pseudobulk
#   3) MEF vs PEA bulk            : sample-level "pseudobulk"
# Output: beta tables + common heatmap (blue-white-red)
# ============================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(decoupleR)
  library(progeny)
  library(readxl)
  library(openxlsx)
  library(pheatmap)
})

# ---- conflict回避（最低限）----
if (requireNamespace("conflicted", quietly=TRUE)) {
  conflicted::conflicts_prefer(dplyr::filter)
  conflicted::conflicts_prefer(dplyr::select)
  conflicted::conflicts_prefer(dplyr::mutate)
  conflicted::conflicts_prefer(dplyr::summarise)
  conflicted::conflicts_prefer(dplyr::arrange)
  conflicted::conflicts_prefer(dplyr::left_join)
}

# -------------------------
# 0) settings
# -------------------------
out_dir <- "C:/Users/myama/Desktop/PROGENy_LRP1assoc_3datasets"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# heatmap color (青-白-赤): あなたの2枚目っぽい見た目に寄せる
lim <- 2
bk  <- seq(-lim, lim, length.out = 101)
cols <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# -------------------------
# 1) utilities
# -------------------------

# (A) pseudobulk mean expression (genes x cells) -> (genes x pseudobulk)
.pseudobulk_mean <- function(expr, group_vec){
  stopifnot(ncol(expr) == length(group_vec))
  grp_levels <- unique(group_vec)
  mm <- Matrix::sparse.model.matrix(~0 + factor(group_vec, levels = grp_levels))
  colnames(mm) <- grp_levels
  n_cells <- Matrix::colSums(mm)
  pb_sum  <- expr %*% mm
  pb_mean <- t(t(pb_sum) / as.numeric(n_cells))
  pb_mean
}

# (B) add PROGENy assay to Seurat (pseudobulk Seurat)
.add_progeny <- function(seu, organism=c("Human","Mouse"), top=1000){
  organism <- match.arg(organism)
  progeny_fun <- get("progeny", envir = asNamespace("decoupleR"))
  seu <- progeny_fun(
    expr = seu,
    organism = organism,
    top = top,
    perm = 1,
    scale = FALSE,
    return_assay = TRUE
  )
  if (!"progeny" %in% Seurat::Assays(seu)) stop("progeny assay が作られていません。")
  seu
}

# (C) LRP1 gene pick (case-insensitive exact first; fallback contains)
.pick_gene <- function(genes, symbol){
  hit <- genes[toupper(genes) == toupper(symbol)]
  if (length(hit) >= 1) return(hit[1])
  hit2 <- genes[grepl(symbol, genes, ignore.case=TRUE)]
  if (length(hit2) >= 1) return(hit2[1])
  stop("Gene not found: ", symbol)
}

# (D) regression beta: pathway ~ LRP1 (+ covariates optional)
# returns data.frame(pathway, beta_LRP1, p_LRP1, adj_p)
.beta_regress <- function(score_mat, lrp1_vec, cov_df=NULL, dataset_name=""){
  # score_mat: pathways x samples
  stopifnot(all(colnames(score_mat) %in% names(lrp1_vec)))
  df_score <- as.data.frame(t(as.matrix(score_mat)))  # samples x pathways
  df_score$pb_id <- rownames(df_score)
  df <- data.frame(pb_id=df_score$pb_id, LRP1=as.numeric(lrp1_vec[df_score$pb_id]), stringsAsFactors=FALSE)
  
  if (!is.null(cov_df)) {
    cov_df <- cov_df[match(df$pb_id, cov_df$pb_id), , drop=FALSE]
    df <- cbind(df, cov_df[, setdiff(colnames(cov_df), "pb_id"), drop=FALSE])
  }
  df <- cbind(df, df_score[, setdiff(colnames(df_score), "pb_id"), drop=FALSE])
  
  pathways <- setdiff(colnames(df_score), "pb_id")
  res <- lapply(pathways, function(pw){
    # covariatesがあれば式に足す
    if (is.null(cov_df)) {
      f <- as.formula(paste0("`", pw, "` ~ LRP1"))
    } else {
      cov_terms <- setdiff(colnames(cov_df), "pb_id")
      f <- as.formula(paste0("`", pw, "` ~ LRP1 + ", paste(cov_terms, collapse=" + ")))
    }
    fit <- lm(f, data=df)
    sm  <- summary(fit)$coefficients
    data.frame(
      dataset = dataset_name,
      pathway = pw,
      beta_LRP1 = sm["LRP1","Estimate"],
      p_LRP1    = sm["LRP1","Pr(>|t|)"],
      stringsAsFactors=FALSE
    )
  }) |> bind_rows() |>
    mutate(adj_p_LRP1 = p.adjust(p_LRP1, method="BH")) |>
    arrange(adj_p_LRP1)
  
  res
}

# -------------------------
# 2) HLCA (fibro_subset) : donor-level pseudobulk -> PROGENy -> beta(LRP1)
# -------------------------
message("=== HLCA fibro_subset ===")

setwd("C:\\Users\\myama\\Desktop\\IP303例")
data_IP <- readRDS("C:\\Users\\myama\\Desktop\\IP303例\\HC_NSIP_ILD_HP_Myo_SSc_303samples.rds")

# biomaRtを使用してENSEMBL IDを遺伝子シンボルに変換
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- rownames(GetAssayData(data_IP, slot = "counts"))
gene_conversion <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl_ids, mart = ensembl)
gene_conversion <- gene_conversion[gene_conversion$hgnc_symbol != "", ]
gene_mapping <- setNames(gene_conversion$hgnc_symbol, gene_conversion$ensembl_gene_id)
counts_data <- GetAssayData(data_IP, slot = "counts")
new_rownames <- ifelse(!is.na(gene_mapping[rownames(counts_data)]), gene_mapping[rownames(counts_data)], rownames(counts_data))
unique_new_rownames <- make.unique(new_rownames)
rownames(counts_data) <- unique_new_rownames
data_data <- GetAssayData(data_IP, slot = "data")
rownames(data_data) <- unique_new_rownames
assay <- CreateAssayObject(counts = counts_data)
assay <- SetAssayData(assay, slot = "data", new.data = data_data)
data_IP[['RNA']] <- assay
head(rownames(GetAssayData(data_IP, slot = "counts")))

# === 線維芽細胞サブセットの抽出 ========================================
fibro_targets <- c("Adventitial fibroblasts", "Alveolar fibroblasts", 
                   "Fibroblasts PLIN2+", "Fibroblasts WIF1+ CHRM2+", 
                   "Peribronchial fibroblasts", "Nerve-associated fibroblasts",
                   "Subpleural fibroblasts", "Secondary crest myofibroblast",
                   "Activated myofibroblasts", "Lipofibroblasts")

fibro_subset <- subset(data_IP, subset = original_ann_level_4 %in% fibro_targets)

stopifnot(exists("fibro_subset"))
hlca <- fibro_subset

# donor id
if (!"donor_id" %in% colnames(hlca@meta.data)) stop("HLCA: donor_id が meta.data にありません。")
donor <- as.character(hlca@meta.data$donor_id)

DefaultAssay(hlca) <- "RNA"
if (ncol(GetAssayData(hlca, slot="data")) == 0) hlca <- NormalizeData(hlca)
expr_hlca <- GetAssayData(hlca, slot="data")  # genes x cells (sparse)

# pseudobulk mean by donor
pb_expr_hlca <- .pseudobulk_mean(expr_hlca, donor)

pb_hlca <- CreateSeuratObject(counts = pb_expr_hlca)
DefaultAssay(pb_hlca) <- "RNA"
pb_hlca <- SetAssayData(pb_hlca, slot="data", new.data = pb_expr_hlca)
pb_hlca@meta.data$pb_id <- colnames(pb_hlca)

# LRP1 vector (per donor pseudobulk)
lrp1_gene_h <- .pick_gene(rownames(pb_expr_hlca), "LRP1")
lrp1_hlca <- as.numeric(pb_expr_hlca[lrp1_gene_h, ])
names(lrp1_hlca) <- colnames(pb_expr_hlca)

# PROGENy scores
pb_hlca <- .add_progeny(pb_hlca, organism="Human", top=1000)
score_hlca <- GetAssayData(pb_hlca, assay="progeny", slot="data") # pathways x donorPB

# beta (HLCA): pathway ~ LRP1
conflicts_prefer(GSEABase::setdiff)
conflicts_prefer(lubridate::setdiff)
conflicts_prefer(GenomicRanges::setdiff)
conflicts_prefer(base::setdiff)
res_hlca <- .beta_regress(score_hlca, lrp1_hlca, cov_df=NULL,
                          dataset_name="HLCA_fibro_pseudobulk")
write.csv(res_hlca, file.path(out_dir, "HLCA_beta_LRP1.csv"), row.names=FALSE)

# -------------------------
# 3) Mouse BLM fibro (makoto.figure3) : mouse-level pseudobulk -> PROGENy -> beta(Lrp1 + day optional)
# -------------------------
message("=== Mouse BLM fibro (makoto.figure3) ===")
setwd("C:\\Users\\myama\\Desktop\\BLM_Seurat")
file_path <- "C:\\Users\\myama\\Desktop\\BLM_Seurat\\BLM_Seurat.qs"  
data_BLM <- qread(file_path)

DimPlot(data_BLM, group.by = "condition") 
DimPlot(data_BLM, group.by = "orig.ident")  
DimPlot(data_BLM, group.by = "celltype") 
DimPlot(data_BLM, group.by = "celltype_init") 
DimPlot(data_BLM, group.by = "celltype_fine") 
DimPlot(data_BLM, group.by = "celltype_coarse") 

DimPlot(data_BLM, group.by = "celltype", label = TRUE, repel = TRUE) + 
  theme_minimal() + 
  theme(text = element_text(size = 1), 
        legend.position = "none",
        panel.grid.major = element_blank(),  # 大きなグリッドを消す
        panel.grid.minor = element_blank())  # 小さなグリッドを消す

#makoto.figure2用UMAP作成
meta <- data_BLM@meta.data
meta$makoto.figure3 <- meta$celltype

## macrophage / monocyte / DC のみ上書きルールを定義
meta$makoto.figure3 <- case_when(
  grepl("^Mac_AM",  meta$celltype_fine) ~ "Mac_Alveolar",      # 肺胞マクロファージ
  grepl("^Mac_int", meta$celltype_fine) |
    grepl("^Mac_IM", meta$celltype_fine) ~ "Mac_Interstitial", # 間質マクロファージ
  grepl("^Mac_proliferated", meta$celltype_fine) ~ "Mac_Interstitial",  # 好みで
  grepl("^Mo_Ly6Chi", meta$celltype_fine) ~ "Monocyte",        # classical 相当
  grepl("^Mo_Ly6Clo", meta$celltype_fine) ~ "Monocyte",        # non-classical 相当
  grepl("^Mo_",       meta$celltype_fine) ~ "Monocyte",        # その他 Mo_ もまとめる
  grepl("^DC_", meta$celltype) | 
    grepl("^DC_", meta$celltype_fine) ~ "DC",
  TRUE ~ meta$celltype)

data_BLM$makoto.figure3 <- meta$makoto.figure3
Idents(data_BLM) <- "makoto.figure3"

cell_order <- c(
  "Epi_AT2", "Epi_AT1", "Epi_ciliated_clara",
  "Fibroblast", "Fibroblast_Peribronchiolar",
  "Endo_capillary", "Endo_arterial", "Endo_venous", "Endo_LEC",
  "Pericyte", "SMC",
  "Mac_Alveolar", "Mac_Interstitial", "Monocyte", "DC",
  "NK", "ILC2",
  "Tcell_CD4T", "Tcell_CD8T", "Tcell_Tgd", "Tcell_Treg",
  "Neutrophil", "Eosinophil", "Basophil",
  "Platelets", "Plasma cell", "Pro-platelet",
  "Sox10pos_cells", "Mesothelial",
  "misc", "doublet"
)

n_clust  <- length(unique(data_BLM$makoto.figure3))
base_cols <- brewer.pal(8, "Set2")
cols <- colorRampPalette(base_cols)(n_clust)
names(cols) <- sort(unique(data_BLM$makoto.figure3))

# stopifnot(exists("data_BLM"))

# fibro抽出: makoto.figure3 の "Fibroblast" & "Fibroblast_Peribronchiolar"
if (!"makoto.figure3" %in% colnames(data_BLM@meta.data)) stop("BLM: makoto.figure3 が meta.data にありません。")
cells_fib_m <- rownames(data_BLM@meta.data)[data_BLM@meta.data$makoto.figure3 %in% c("Fibroblast","Fibroblast_Peribronchiolar")]
if (length(cells_fib_m) == 0) stop("BLM: fibro cells が 0 です。makoto.figure3 のラベルを確認してください。")
blm_fib <- subset(data_BLM, cells = cells_fib_m)

# sample (mouse) id
if (!"orig.ident" %in% colnames(blm_fib@meta.data)) stop("BLM: orig.ident がありません。")
mouse_id <- as.character(blm_fib@meta.data$orig.ident)

# day (optional covariate): condition から数値化（あなたのデータは day00/day10/day28 など）
cond_col <- if ("condition" %in% colnames(blm_fib@meta.data)) "condition" else stop("BLM: condition 列がありません。")
day_num <- suppressWarnings(as.numeric(sub(".*?(\\d+).*", "\\1", as.character(blm_fib@meta.data[[cond_col]]))))
# NAが多いなら covariate入れないほうが安全
use_day_cov <- sum(!is.na(day_num)) >= 0.8*nrow(blm_fib@meta.data)

DefaultAssay(blm_fib) <- "RNA"
if (ncol(GetAssayData(blm_fib, slot="data")) == 0) blm_fib <- NormalizeData(blm_fib)
expr_blm <- GetAssayData(blm_fib, slot="data")

# pseudobulk mean by mouse_id
pb_expr_blm <- .pseudobulk_mean(expr_blm, mouse_id)

pb_blm <- CreateSeuratObject(counts = pb_expr_blm)
DefaultAssay(pb_blm) <- "RNA"
pb_blm <- SetAssayData(pb_blm, slot="data", new.data = pb_expr_blm)
pb_blm@meta.data$pb_id <- colnames(pb_blm)

# Lrp1 vector
lrp1_gene_m <- .pick_gene(rownames(pb_expr_blm), "Lrp1")
lrp1_blm <- as.numeric(pb_expr_blm[lrp1_gene_m, ])
names(lrp1_blm) <- colnames(pb_expr_blm)

# covariate day per pseudobulk（mouse単位で平均 day を持たせる：mouse内で複数dayが混ざる場合の保険）
cov_blm <- NULL
if (use_day_cov) {
  md <- blm_fib@meta.data
  md$mouse_id <- mouse_id
  md$day_num  <- day_num
  day_by_mouse <- md %>%
    group_by(mouse_id) %>%
    summarise(day_num = mean(day_num, na.rm=TRUE), .groups="drop")
  cov_blm <- data.frame(pb_id = day_by_mouse$mouse_id, day_num = day_by_mouse$day_num)
}

# PROGENy
pb_blm <- .add_progeny(pb_blm, organism="Mouse", top=1000)
score_blm <- GetAssayData(pb_blm, assay="progeny", slot="data")

# beta (BLM): pathway ~ Lrp1 (+ day_num if usable)
res_blm <- .beta_regress(score_blm, lrp1_blm, cov_df=cov_blm,
                         dataset_name = if (use_day_cov) "MouseBLM_fibro_pseudobulk_dayAdj" else "MouseBLM_fibro_pseudobulk")
write.csv(res_blm, file.path(out_dir, "MouseBLM_beta_Lrp1.csv"), row.names=FALSE)

# -------------------------
# 4) MEF vs PEA bulk : PROGENy -> beta(Lrp1)   (※ここは「サンプルが少ない」ので解釈注意)
# -------------------------
message("=== MEF vs PEA bulk ===")

xlsx_m <- "C:/Users/myama/Desktop/MEF.PEA-13.RNAseq結果.xlsx"

mraw <- readxl::read_excel(xlsx_m, sheet = 1) |> as.data.frame()

# mraw <- readxl::read_excel(xlsx_m, sheet=1) |> as.data.frame()
stopifnot(exists("mraw"))

gene_col <- colnames(mraw)[1]
genes <- as.character(mraw[[gene_col]])
expr_cols <- grep("MEF_RNA|PEA_RNA", colnames(mraw), value=TRUE)
if (length(expr_cols) < 2) stop("MEF/PEA: MEF_RNA/PEA_RNA列が見つかりません。")

expr_mef <- as.matrix(mraw[, expr_cols, drop=FALSE])
rownames(expr_mef) <- genes
expr_mef <- apply(expr_mef, 2, function(x) as.numeric(as.character(x)))
rownames(expr_mef) <- genes
colnames(expr_mef) <- make.unique(colnames(expr_mef))

pb_mef <- CreateSeuratObject(counts = expr_mef)
DefaultAssay(pb_mef) <- "RNA"
pb_mef <- SetAssayData(pb_mef, slot="data", new.data = as.matrix(expr_mef))
pb_mef@meta.data$pb_id <- colnames(pb_mef)

lrp1_gene_mef <- .pick_gene(rownames(expr_mef), "Lrp1")
lrp1_mef <- as.numeric(expr_mef[lrp1_gene_mef, ])
names(lrp1_mef) <- colnames(expr_mef)

pb_mef <- .add_progeny(pb_mef, organism="Mouse", top=1000)
score_mef <- GetAssayData(pb_mef, assay="progeny", slot="data")

res_mef <- .beta_regress(score_mef, lrp1_mef, cov_df=NULL,
                         dataset_name="MEF_PEA_bulk")
write.csv(res_mef, file.path(out_dir, "MEF_PEA_beta_Lrp1.csv"), row.names=FALSE)

# -------------------------
# 5) Combine -> common pathways -> heatmap (z-scored by dataset)
# -------------------------
message("=== Combine -> heatmap ===")

eff_all <- bind_rows(
  res_hlca %>% transmute(dataset, pathway, effect = beta_LRP1),
  res_blm  %>% transmute(dataset, pathway, effect = beta_LRP1),
  res_mef  %>% transmute(dataset, pathway, effect = beta_LRP1)
)

write.csv(eff_all, file.path(out_dir, "ALL_beta_long.csv"), row.names=FALSE)

# common pathways in all 3
conflicts_prefer(matrixStats::count)
conflicts_prefer(dplyr::count)
common_pw <- eff_all %>%
  count(pathway) %>%
  filter(n == n_distinct(eff_all$dataset)) %>%
  pull(pathway)

eff_common <- eff_all %>%
  filter(pathway %in% common_pw)

mat <- eff_common %>%
  pivot_wider(names_from = pathway, values_from = effect) %>%
  as.data.frame()

rownames(mat) <- mat$dataset
mat$dataset <- NULL
mat <- as.matrix(mat)

# z-score by dataset (row-wise)
mat_z <- t(scale(t(mat)))

# clip to [-2,2] for stable colors
mat_z2 <- mat_z
mat_z2[mat_z2 >  lim] <-  lim
mat_z2[mat_z2 < -lim] <- -lim

write.csv(mat,   file.path(out_dir, "HEATMAP_matrix_beta.csv"))
write.csv(mat_z, file.path(out_dir, "HEATMAP_matrix_beta_z.csv"))

pdf(file.path(out_dir, "HEATMAP_common_PROGENy_beta_LRP1.pdf"),
    width = 10, height = max(4, 0.35*nrow(mat_z2)+2))

lim <- 2
bk  <- seq(-lim, lim, length.out = 101)
cols <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

mat_z2 <- mat_z
mat_z2[mat_z2 >  lim] <-  lim
mat_z2[mat_z2 < -lim] <- -lim

pheatmap::pheatmap(
  mat_z2,
  color = cols,
  breaks = bk,          # ← これ超重要（中心=白に固定）
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Common PROGENy pathways associated with LRP1 (beta; z-scored by dataset)"
)

dev.off()

tiff(file.path(out_dir, "HEATMAP_common_PROGENy_beta_LRP1_1200dpi.tiff"),
     width=12, height=5, units="in", res=1200, compression="lzw")
lim <- 2
bk  <- seq(-lim, lim, length.out = 101)
cols <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

mat_z2 <- mat_z
mat_z2[mat_z2 >  lim] <-  lim
mat_z2[mat_z2 < -lim] <- -lim

pheatmap::pheatmap(
  mat_z2,
  color = cols,
  breaks = bk,          # ← これ超重要（中心=白に固定）
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Common PROGENy pathways associated with LRP1 (beta; z-scored by dataset)"
)

dev.off()

message("DONE. outputs in: ", out_dir)


