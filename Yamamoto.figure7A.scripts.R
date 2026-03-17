########Yamamoto.figure7A.scripts.R###########

data_dir <- "data/"
integrated_obj <- readRDS(file.path(data_dir, "seuratObject.RDS"))

# --- 今回使うFOV（WT21-2, C3, W6） ---
target_fovs <- c(9:17, 27:35, 36:43)

cells_to_keep <- rownames(integrated_obj@meta.data[integrated_obj@meta.data$fov %in% target_fovs, ])
obj <- subset(integrated_obj, cells = cells_to_keep)

# --- sample_name / lung_group ---
obj$sample_name <- NA
obj$sample_name[obj$fov %in% 9:17]   <- "WT21-2"
obj$sample_name[obj$fov %in% 27:35]  <- "C3"
obj$sample_name[obj$fov %in% 36:43]  <- "W6"

obj$lung_group <- NA
obj$lung_group[obj$sample_name %in% c("C3","W6")] <- "WT_Unstimulated"
obj$lung_group[obj$sample_name == "WT21-2"]       <- "WT_Stimulated"
obj$lung_group <- factor(obj$lung_group, levels = c("WT_Unstimulated","WT_Stimulated"))

DimPlot(obj, group.by="lung_group", label=TRUE)
VlnPlot(obj, features="LRP1", group.by="lung_group", pt.size=0)

#QC～SCT
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
VlnPlot(obj, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
obj <- subset(obj, subset = nFeature_RNA > 100 & nFeature_RNA < 600)
obj <- SCTransform(obj, verbose = TRUE)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:30)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
DimPlot(obj, reduction = "umap", label = TRUE)

#H&Eで見た情報をmeta.dataに追加
md <- obj@meta.data
md$fov <- as.integer(md$fov)
md$region <- NA

# BLM側
md$region[md$fov %in% 9:11]  <- "BLM_central_highFib"
md$region[md$fov %in% 12:14] <- "BLM_mid_midFib"
md$region[md$fov %in% 15:17] <- "BLM_distal_subpleural_mixed"

# 正常 W6
md$region[md$fov %in% c(36,37,38,39,42,43)] <- "Ctrl_bronchial_wall"
md$region[md$fov %in% c(40,41)]             <- "Ctrl_distal_alveolar"

# 正常 C3
md$region[md$fov %in% c(27,28,31,33)] <- "Ctrl_central"
md$region[md$fov %in% c(29,32,34)]    <- "Ctrl_mid"
md$region[md$fov %in% c(30,35)]       <- "Ctrl_distal_subpleural"

obj@meta.data <- md
VlnPlot(obj, features="LRP1", group.by="region", pt.size=0) + NoLegend()

#SinglRでアノテーション
lung_h5ad <- readH5AD("C://Users//yamas//OneDrive//デスクトップ//RStudioの資料そのほか//tabula-muris-senis-droplet-processed-official-annotations-Lung.h5ad")
test_sce <- as.SingleCellExperiment(obj)
sct_logdata <- GetAssayData(obj, slot = "data", assay = "SCT")
logcounts(test_sce) <- as.matrix(sct_logdata)
logcounts(lung_h5ad) <- assay(lung_h5ad, "X")

singleR_res <- SingleR(
  test = test_sce,
  ref = lung_h5ad,
  labels = lung_h5ad$free_annotation,
  assay.type.test = "logcounts"
)

obj$SingleR_label <- singleR_res$labels
DimPlot(obj, group.by = "SingleR_label", label = TRUE, repel = TRUE)

#SCTデータで線維化スコア
DefaultAssay(obj) <- "SCT"
fib_genes <- c("Col1a1","Col1a2","Col3a1","Fn1","Acta2","Tagln","Ctgf","Postn","Lox")
# 1) 存在する遺伝子だけに絞る
conflicts_prefer(Biostrings::intersect)
conflicts_prefer(GSEABase::intersect)
conflicts_prefer(lubridate::intersect)
conflicts_prefer(base::intersect)
present <- intersect(fib_genes, rownames(obj[["SCT"]]))
present
length(present)

# 2) AddModuleScore（ctrlを小さくする）
obj <- AddModuleScore(
  obj,
  features = list(present),
  name = "FibScore",
  ctrl = 20   # ← まず20。まだ落ちるなら 10, 5 と下げる
)

find_seg_image <- function(fov_id, type = "Cell"){
  base <- "E:/Slide1_20250305_18_03_2025_6_39_17_643"
  tifs <- list.files(base, recursive=TRUE, pattern="(?i)\\.tif(f)?$", full.names=TRUE)
  
  pattern <- if(type == "Cell"){
    sprintf("CellLabels_F%05d", fov_id)
  } else {
    sprintf("CompartmentLabels_F%05d", fov_id)
  }
  
  cand <- tifs[grepl(pattern, tifs)]
  cand
}

#画像の構成
plot_overlay_seg <- function(obj, fov_id, type="Cell"){
  
  img_list <- find_seg_image(fov_id, type)
  if(length(img_list)==0){
    stop("Segmentation image not found for FOV ", fov_id)
  }
  
  img_path <- img_list[1]
  message("Using image: ", img_path)
  
  obj_fov <- subset(obj, subset = fov == fov_id)
  md <- obj_fov@meta.data
  
  # 発現データをmeta.dataに追加（←ここだけ LRP1 に）
  expr <- FetchData(obj_fov, vars = "LRP1")
  md$LRP1 <- expr$LRP1
  
  img <- tiff::readTIFF(img_path, native=TRUE)
  g <- grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
  
  ggplot(md, aes(x = x_FOV_px, y = y_FOV_px)) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(aes(color = LRP1), size=0.3, alpha=0.7) +
    scale_color_viridis_c(option="C") +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("FOV ", fov_id, " – LRP1 (", type, " labels)"))
}


p9  <- plot_overlay_seg(obj, 9)
p13 <- plot_overlay_seg(obj, 13)
p16 <- plot_overlay_seg(obj, 16)
p9 + p13 + p16

#FOVごとにLRP1の発現を定量化
DefaultAssay(obj) <- "SCT"  # 念のため（SCTで見たいなら）
obj$LRP1_expr <- FetchData(obj, vars = "LRP1")[,1]

df_fov <- obj@meta.data %>%
  group_by(fov) %>%
  summarise(LRP1_mean = mean(LRP1_expr, na.rm = TRUE),
            LRP1_pos  = mean(LRP1_expr > 0, na.rm = TRUE) * 100,
            n_cells   = dplyr::n(),
            .groups = "drop")

df_fov

obj$WT21_region <- NA
obj$WT21_region[obj$fov %in% 9:11]  <- "Central (severe)"
obj$WT21_region[obj$fov %in% 12:14] <- "Middle"
obj$WT21_region[obj$fov %in% 15:17] <- "Peripheral (subpleural)"

md_wt21 <- dplyr::filter(obj@meta.data, fov %in% 9:17)
df_reg <- md_wt21 %>%
  group_by(WT21_region) %>%
  summarise(LRP1_mean = mean(LRP1_expr, na.rm=TRUE),
            LRP1_pos  = mean(LRP1_expr > 0, na.rm=TRUE) * 100,
            .groups="drop")

df_reg

ggplot(df_reg, aes(x = WT21_region, y = LRP1_pos)) +
  geom_col() +
  theme_classic() +
  xlab("") + ylab("LRP1+ cells (%)") +
  theme(axis.text.x = element_text(angle=20, hjust=1))


obj$WT21_pathology <- NA
obj$WT21_pathology[obj$fov %in% 9:14] <- "Fibrotic core (FOV9–14)"
obj$WT21_pathology[obj$fov %in% 15:17] <- "Peripheral / subpleural (FOV15–17)"

md_wt21b <- dplyr::filter(obj@meta.data, fov %in% 9:17)

df_core <- md_wt21b %>%
  dplyr::group_by(WT21_pathology) %>%
  dplyr::summarise(
    LRP1_mean = mean(LRP1_expr, na.rm=TRUE),
    LRP1_pos  = mean(LRP1_expr > 0, na.rm=TRUE) * 100,
    .groups="drop"
  )
df_core

ggplot(df_core, aes(x = WT21_pathology, y = LRP1_pos)) +
  geom_col() +
  theme_classic() +
  xlab("") + ylab("LRP1+ cells (%)") +
  theme(axis.text.x = element_text(angle=15, hjust=1))

md_wt21_all <- dplyr::filter(obj@meta.data, fov %in% 9:17)
vmax <- as.numeric(quantile(md_wt21_all$LRP1_expr, 0.99, na.rm=TRUE))

scale_color_viridis_c(option="C", limits=c(0, vmax), oob=scales::squish)

plot_overlay_seg_fixed <- function(obj, fov_id, vmax, type="Cell"){
  
  img_list <- find_seg_image(fov_id, type)
  if(length(img_list)==0) stop("Segmentation image not found for FOV ", fov_id)
  
  img_path <- img_list[1]
  message("Using image: ", img_path)
  
  obj_fov <- subset(obj, subset = fov == fov_id)
  md <- obj_fov@meta.data
  
  # LRP1をmeta.dataに引っ張る（SCT/RNAどちらでもOK）
  expr <- FetchData(obj_fov, vars = "LRP1")[,1]
  md$LRP1_plot <- expr
  
  img <- tiff::readTIFF(img_path, native=TRUE)
  g <- grid::rasterGrob(img, width=unit(1,"npc"), height=unit(1,"npc"))
  
  ggplot(md, aes(x = x_FOV_px, y = y_FOV_px)) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(aes(color = LRP1_plot), size=1.5, alpha=0.7) +
    scale_color_viridis_c(option="C", limits=c(0, vmax), oob=scales::squish) +
    coord_fixed() +
    scale_y_reverse() +
    theme_void() +
    ggtitle(paste0("FOV ", fov_id, " – LRP1"))
}

p9  <- plot_overlay_seg_fixed(obj, 9,  vmax)
p13 <- plot_overlay_seg_fixed(obj, 13, vmax)
p16 <- plot_overlay_seg_fixed(obj, 16, vmax)

library(patchwork)
p_overlay <- p9 + p13 + p16 + plot_layout(guides = "collect") & theme(legend.position="right")
p_overlay

outfile <- file.path(save_dir, "LRP1_overlay_FOV9_13_16.tiff")

ggsave(
  filename = outfile,
  plot = p_overlay,
  device = "tiff",
  width = 10,
  height = 4,
  units = "in",
  dpi = 1200,
  compression = "lzw"
)


obj$WT_group <- NA
obj$WT_group[obj$fov %in% 9:17] <- "WT21 (BLM d21)"
obj$WT_group[obj$fov %in% c(27:35,36:43)] <- "WT Unstim"

md_wt_compare <- dplyr::filter(obj@meta.data, WT_group %in% c("WT21 (BLM d21)","WT Unstim"))

df_wtcomp <- md_wt_compare %>%
  dplyr::group_by(WT_group) %>%
  dplyr::summarise(LRP1_pos = mean(LRP1_expr > 0)*100,
                   LRP1_mean = mean(LRP1_expr),
                   .groups="drop")

df_wtcomp

# WT21: 9–14 を線維化コア、15–17 を末梢として定義
obj$WT21_pathology <- NA
obj$WT21_pathology[obj$fov %in% 9:14] <- "Fibrotic core (FOV9-14)"
obj$WT21_pathology[obj$fov %in% 15:17] <- "Peripheral/subpleural (FOV15-17)"

md_wt21b <- dplyr::filter(obj@meta.data, fov %in% 9:17)

df_core <- md_wt21b %>%
  dplyr::group_by(WT21_pathology) %>%
  dplyr::summarise(
    LRP1_pos  = mean(LRP1_expr > 0, na.rm=TRUE) * 100,
    LRP1_mean = mean(LRP1_expr, na.rm=TRUE),
    n_cells   = dplyr::n(),
    .groups="drop"
  )

df_core

#保存
# ===== 保存先ディレクトリの指定 =====
save_dir <- "C:/Users/yamas/OneDrive/デスクトップ/CosMxマウス解析"

# フォルダが無ければ作成
if(!dir.exists(save_dir)){
  dir.create(save_dir, recursive = TRUE)
}

# ===== CSVの保存 =====
write.csv(df_core, file.path(save_dir, "WT21_LRP1_core_vs_peripheral.csv"), row.names = FALSE)
write.csv(df_wtcomp, file.path(save_dir, "WT_LRP1_unstim_vs_BLMd21.csv"), row.names = FALSE)
write.csv(df_fov, file.path(save_dir, "WT_LRP1_byFOV.csv"), row.names = FALSE)
# ===== long形式データの保存 =====
md_wt21b_long <- md_wt21b %>%
  dplyr::filter(!is.na(WT21_pathology)) %>%
  dplyr::select(fov, WT21_pathology, LRP1_expr)

write.csv(md_wt21b_long, file.path(save_dir, "WT21_cells_LRP1expr_core_vs_peripheral_long.csv"), row.names = FALSE)

md_wt_compare_long <- md_wt_compare %>%
  dplyr::select(fov, WT_group, LRP1_expr)

write.csv(md_wt_compare_long, file.path(save_dir, "WT_cells_LRP1expr_unstim_vs_BLMd21_long.csv"), row.names = FALSE)

#Fibrosis scoreを調べる
grep("^FibScore", colnames(obj@meta.data), value = TRUE)
fib_col <- grep("^FibScore", colnames(obj@meta.data), value=TRUE)[1]

df_fib_fov <- obj@meta.data %>%
  dplyr::group_by(fov) %>%
  dplyr::summarise(
    FibScore_mean = mean(.data[[fib_col]], na.rm=TRUE),
    FibScore_med  = median(.data[[fib_col]], na.rm=TRUE),
    n_cells       = dplyr::n(),
    .groups="drop"
  )
df_fib_fov

ggplot(df_fib_fov, aes(x = factor(fov), y = FibScore_mean)) +
  geom_col() +
  theme_classic() +
  xlab("FOV") + ylab("Fibrosis module score (mean)") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

fib_col <- grep("^FibScore", colnames(obj@meta.data), value=TRUE)[1]

md_unstim <- dplyr::filter(obj@meta.data, fov %in% c(27:35, 36:43))
thr_fib <- as.numeric(quantile(md_unstim[[fib_col]], 0.95, na.rm=TRUE))
thr_fib


df_fibpct_fov <- obj@meta.data %>%
  dplyr::mutate(FibroticCell = .data[[fib_col]] >= thr_fib) %>%
  dplyr::group_by(fov) %>%
  dplyr::summarise(
    Fibrotic_pct = mean(FibroticCell, na.rm=TRUE) * 100,
    n_cells = dplyr::n(),
    .groups="drop"
  )
df_fibpct_fov

# すでにある df_fov (LRP1_mean/LRP1_pos) と結合
df_merge <- dplyr::left_join(df_fov, df_fib_fov, by="fov")

ggplot(df_merge, aes(x = FibScore_mean, y = LRP1_pos)) +
  geom_point(size=3) +
  theme_classic() +
  xlab("Fibrosis score (mean, per FOV)") +
  ylab("LRP1+ cells (%)")

write.csv(df_fib_fov, file.path(save_dir, "WT_FibScore_byFOV.csv"), row.names=FALSE)
write.csv(df_fibpct_fov, file.path(save_dir, "WT_FibroticPct_byFOV.csv"), row.names=FALSE)

df_merge <- dplyr::left_join(df_fov, df_fib_fov, by="fov")
write.csv(df_merge, file.path(save_dir, "WT_FOV_Fibrosis_vs_LRP1_merged.csv"), row.names=FALSE)
