#Figure3A

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

p_fig3 <- DimPlot(
  data_BLM,
  group.by = "makoto.figure3",
  order = cell_order,
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  scale_color_manual(values = cols) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    
    # ↓ ここで軸の線・目盛り・ラベルを全部消す
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "BLM lung scRNAseq (makoto.figure3)")

p_fig3

ggsave(
  filename = "BLM_lung_scRNAseq_makoto_figure3.tiff",
  plot     = p_fig3,
  device   = "tiff",
  dpi      = 1200,      # 解像度
  width    = 18,        # 図の横幅（cm 単位でも inch でも可）
  height   = 9,         # 図の縦幅
  units    = "cm",
  compression = "lzw",
  limitsize   = FALSE   # 大きな図でもエラーにしない
)


# LRP1の遺伝子発現データのプロット グレースケール
FeaturePlot(data_BLM, features = "Lrp1", min.cutoff = "q10", max.cutoff = "q60", cols = c("lightgrey", "black"), pt.size = 0.2) +
  scale_color_gradientn(colors = c("lightgrey", "darkgrey", "black")) +  # グレースケールに設定
  theme_minimal() +
  labs(title = "Lrp1 Expression in Fibroblasts",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  # 格子を消す


unique_celltypes <- unique(data_BLM@meta.data$celltype)
print(unique_celltypes)

unique_celltypes <- unique(data_BLM@meta.data$orig.ident)
print(unique_celltypes)


#Figure3B
#2025年12月dot plot violin plotの作成
# day00 → Control、それ以外 → BLM
data_BLM$BLM_group <- ifelse(data_BLM$day_group == "day00", "Control", "BLM")
cell_order <- c(
  "Bcell", "DC", "Endo_capillary",
  "Epi_AT1", "Epi_AT2", "Epi_ciliated_clara",
  "Fibroblast", "Fibroblast_Peribronchiolar",
  "MoMac", "Mac_Alveolar", "Mac_Interstitial", "Monocyte", "Pericyte", "SMC"
)

filtered_cells <- data_BLM$makoto.figure3 %in% cell_order
filtered_data_BLM <- subset(data_BLM, subset = makoto.figure3 %in% cell_order)

filtered_data_BLM$BLM_group <- factor(filtered_data_BLM$BLM_group,
                                      levels = c("Control", "BLM"))

VlnPlot(filtered_data_BLM,
        features = "Lrp1",
        group.by = "makoto.figure3",
        split.by = "BLM_group",
        split.plot = TRUE,
        pt.size = 0) +
  scale_fill_manual(values = c("Control" = "grey80", "BLM" = "grey30")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "LRP1 Expression Across Cell Types (Control Left)")

tiff("LRP1_violin_plot.tiff",
     width = 3000, height = 2000, res = 600, compression = "lzw")

# ここに作図コード
VlnPlot(filtered_data_BLM,
        features = "Lrp1",
        group.by = "makoto.figure3",
        split.by = "BLM_group",
        split.plot = TRUE,
        pt.size = 0) +
  scale_fill_manual(values = c("Control" = "grey80", "BLM" = "grey30")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

dev.off()



DotPlot(
  filtered_data_BLM,
  features = "Lrp1",
  group.by = "makoto.figure3",
  cols = c("white", "black"),
  dot.scale = 10
) +
  ggtitle("LRP1 Expression Across Selected Cell Types") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

tiff("LRP1_DotPlot_highres.tiff",
     width = 3500, height = 3000,
     res = 1200, compression = "lzw")

DotPlot(
  filtered_data_BLM,
  features = "Lrp1",
  group.by = "makoto.figure3",
  cols = c("white", "black"),
  dot.scale = 10
) +
  ggtitle("LRP1 Expression Across Selected Cell Types") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

dev.off()


#LRP1の統計検定
lrp1_name <- "Lrp1"
celltypes <- unique(filtered_data_BLM$makoto.figure3)
p_values <- data.frame(celltype = character(), p_value = numeric())

for (ct in celltypes) {
  
  control_cells <- WhichCells(
    filtered_data_BLM,
    expression = (BLM_group == "Control" & makoto.figure3 == ct)
  )
  
  blm_cells <- WhichCells(
    filtered_data_BLM,
    expression = (BLM_group == "BLM" & makoto.figure3 == ct)
  )
  
  if (length(control_cells) > 0 & length(blm_cells) > 0) {
    
    control_data <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = control_cells)[, 1]
    blm_data     <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = blm_cells)[, 1]
    
    p_val <- wilcox.test(control_data, blm_data)$p.value
    
    p_values <- rbind(
      p_values,
      data.frame(celltype = ct, p_value = p_val)
    )
  }
}

p_values

#LRP1が優位に高いクラスターは
lrp1_name <- "Lrp1"

increase_list <- data.frame(celltype = character(), 
                            mean_control = numeric(),
                            mean_BLM = numeric(),
                            stringsAsFactors = FALSE)

for (ct in p_values$celltype[p_values$p_value < 0.05]) {
  
  control_cells <- WhichCells(filtered_data_BLM, expression = (BLM_group == "Control" & makoto.figure3 == ct))
  blm_cells     <- WhichCells(filtered_data_BLM, expression = (BLM_group == "BLM" & makoto.figure3 == ct))
  
  control_mean <- mean(FetchData(filtered_data_BLM, vars = lrp1_name, cells = control_cells)[[lrp1_name]])
  blm_mean     <- mean(FetchData(filtered_data_BLM, vars = lrp1_name, cells = blm_cells)[[lrp1_name]])
  
  if (blm_mean > control_mean) {
    increase_list <- rbind(increase_list, 
                           data.frame(celltype = ct, 
                                      mean_control = control_mean, 
                                      mean_BLM = blm_mean))
  }
}

increase_list


#Figure3C
#コネクトーム解析

unique(data_BLM@meta.data$condition)

# グループの対応表
condition_grouping <- list(
  "day00UT" = "day00",
  "day03BLM3mg" = "day03",
  "day03BLM1.25mg" = "day03",
  "day07BLM1.25mg" = "day07",
  "day07BLM3mg" = "day07",
  "day14BLM1.25mg" = "day14",
  "day14BLM3mg" = "day14",
  "day28BLM1.25mg" = "day28",
  "day42BLM1.25mg" = "day42"
)

# Seuratオブジェクトに新しい列 "condition_grouped" を追加
data_BLM@meta.data$condition_grouped <- plyr::mapvalues(
  data_BLM@meta.data$condition,
  from = names(condition_grouping),
  to   = base::unname(condition_grouping)
)

# day63 を除く
data_filtered <- subset(data_BLM, subset = condition_grouped != "day63")

# 分割
data_by_day <- SplitObject(data_filtered, split.by = "condition_grouped")

# 出力を格納するリスト
connectome_by_day <- list()

# 各日付でループ
conflicted::conflicts_prefer(scales::rescale)

for (day in names(data_by_day)) {
  message("Processing: ", day)
  
  # 細胞分類を 'celltype_coarse' に設定
  Idents(data_by_day[[day]]) <- data_by_day[[day]]@meta.data$celltype_coarse
  
  # 正規化とスケーリング（未処理であれば）
  data_by_day[[day]] <- NormalizeData(data_by_day[[day]])
  data_by_day[[day]] <- ScaleData(data_by_day[[day]])
  
  # Connectome の作成
  connectome_by_day[[day]] <- CreateConnectome(data_by_day[[day]], species = "mouse", p.values = TRUE)
}


# 可視化したい細胞群
cells_of_interest <- c("Mac_AM", "Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar")

# 各日付で CircosPlot を作成
conflicted::conflicts_prefer(clusterProfiler::select)
conflicted::conflicts_prefer(ComplexHeatmap::draw)


for (day in names(connectome_by_day)) {
  message("Plotting Circos: ", day)
  
  CircosPlot(connectome_by_day[[day]],
             sources.include = cells_of_interest,
             targets.include = cells_of_interest,
             min.pct = 0.1,
             min.z = 1,
             lab.cex = 0.6,
             title = paste("Connectome Circos Plot -", day))
}

connectome_list <- connectome_by_day


# ソース（リガンド放出側）を Mac_AM のみに限定
source_cell <- "Mac_AM"

# ターゲット（リガンド受容側）の細胞群
target_cells <- c("Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar")

# 各日付で CircosPlot を作成
conflicted::conflicts_prefer(clusterProfiler::select)
conflicted::conflicts_prefer(ComplexHeatmap::draw)

for (day in names(connectome_by_day)) {
  message("Plotting Circos: ", day)
  
  CircosPlot(connectome_by_day[[day]],
             sources.include = source_cell,          # Mac_AM のみ
             targets.include = target_cells,         # Fibro 系の細胞群
             min.pct = 0.1,
             min.z = 1,
             lab.cex = 0.6,
             title = paste("Connectome Circos Plot -", day))
}

connectome_list <- connectome_by_day


for (day in names(connectome_by_day)) {
  message("Plotting Circos: ", day)
  
  CircosPlot(connectome_by_day[[day]],
             sources.include = source_cell,
             targets.include = target_cells,
             min.pct = 0.1,
             min.z = 1,
             lab.cex = 1.5,  # 文字の大きさ
             title = paste("Connectome Circos Plot -", day))
}

connectome_list <- connectome_by_day

# Tgfb1を用いた長期的相互作用の可視化
EdgePlot <- function(connectome.list, 
                     mech = "Tgfb1", 
                     sources.use = NULL, 
                     targets.use = NULL) {
  
  require(dplyr)
  require(ggplot2)
  
  combined <- connectome.list[[1]] %>%
    dplyr::filter(ligand == mech)
  combined$sample <- names(connectome.list)[1]
  
  for (i in 2:length(connectome.list)) {
    tmp <- connectome.list[[i]] %>%
      dplyr::filter(ligand == mech)
    tmp$sample <- names(connectome.list)[i]
    combined <- rbind(combined, tmp)
  }
  
  if (!is.null(sources.use)) {
    combined <- combined %>% dplyr::filter(source %in% sources.use)
  }
  if (!is.null(targets.use)) {
    combined <- combined %>% dplyr::filter(target %in% targets.use)
  }
  
  p <- ggplot(combined, aes(x = sample, y = weight_norm, group = target, color = target)) +
    geom_point(size = 2) +
    geom_line() +
    facet_wrap(~source, scales = "free_y") +
    theme_minimal() +
    labs(title = paste0("Longitudinal Interaction via ", mech),
         y = "Normalized Interaction Weight",
         x = "Timepoint")
  
  return(p)
}


Longitudinal <- EdgePlot(
  connectome_list,
  mech = "Tgfb1",
  sources.use = "Mac_AM",
  targets.use = "Fibro_adven"
)

print(Longitudinal)



# " - Lrp1" または "Lrp1 - " を含むペアをフィルタリング
lrp1_related_pairs <- pair_list[grepl("Lrp1", pair_list)]

# 確認（先頭10個など）
head(sort(lrp1_related_pairs), 10)


# connectome_list から必要な情報を抽出
conflicted::conflicts_prefer(stats::cor)

library(dplyr)

tgfb1_df <- do.call(rbind, lapply(names(connectome_list), function(day) {
  df <- connectome_list[[day]]
  df$Timepoint <- day
  return(df)
})) %>%
  dplyr::filter(ligand == "Tgfb1", source == "Mac_AM", target == "Fibro_adven") %>%
  dplyr::mutate(pair = paste0(ligand, "-", receptor))

library(ggplot2)

ggplot(tgfb1_df, aes(x = Timepoint, y = weight_norm, fill = pair)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Tgfb1-mediated signaling: Mac_AM → Fibro_adven",
       x = "Timepoint", y = "Normalized Interaction Weight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Lrp1を含む行だけを抽出したConnectomeデータを日付ごとに準備
connectome_lrp1_only <- lapply(connectome_by_day, function(df) {
  df[grepl("Lrp1", df$pair), ]
})

# 可視化したい細胞群
cells_of_interest <- c("Mac_AM", "Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar")

# Lrp1限定 CircosPlot を日付ごとに作成
for (day in names(connectome_lrp1_only)) {
  message("Plotting Lrp1 Circos: ", day)
  
  CircosPlot(connectome_lrp1_only[[day]],
             sources.include = cells_of_interest,
             targets.include = cells_of_interest,
             min.pct = 0.05,     # 必要なら閾値調整
             min.z = 0.5,
             lab.cex = 1,
             title = paste("Lrp1-related Connectome -", day))
}



# day42 の Connectome データを取得
df_day42 <- connectome_by_day[["day42"]]

# Lrp1 が受容体、かつ target が Fibro_peribron の行のみ抽出
lrp1_fibro_peribron <- df_day42[df_day42$receptor == "Lrp1" & df_day42$target == "Fibro_peribron", ]

# 結果を確認
head(lrp1_fibro_peribron)

# 可視化のために抽出したデータだけを使って CircosPlot を作成
CircosPlot(
  lrp1_fibro_peribron,
  sources.include = unique(lrp1_fibro_peribron$source),  # 自動でソースを設定
  targets.include = "Fibro_peribron",
  min.pct = 0.05,   # 必要に応じて調整
  min.z = 0.5,
  lab.cex = 1.2,
  title = "Lrp1 signaling to Fibro_peribron (day42)"
)

unique(data_BLM@meta.data$celltype_coarse)

#２０２５年７月１４日
# マクロファージと線維芽細胞のクラスをまとめ、それ以外は既存のラベルを保持
data_BLM@meta.data$celltype_coarse_AMfibro_merged <- dplyr::case_when(
  data_BLM@meta.data$celltype_coarse %in% c("Mac_AM", "Mac_IM", "Mac_int", "Mac_proliferated") ~ "Mac_AM_group",
  data_BLM@meta.data$celltype_coarse %in% c("Fibro_adven", "Fibro_alveolar", "Fibro_peribron", "Fibro_proliferated") ~ "Fibroblast_group",
  TRUE ~ data_BLM@meta.data$celltype_coarse  # それ以外は元のまま
)

DimPlot(data_BLM, group.by = "celltype_coarse_AMfibro_merged", label = TRUE, repel = TRUE) +
  ggtitle("Merged Macrophage and Fibroblast Classes")






original_meta <- data_BLM@meta.data[, "celltype_coarse_AMfibro_merged", drop = FALSE]
split_objs <- SplitObject(data_BLM, split.by = "condition_grouped")

for (i in names(split_objs)) {
  cells <- colnames(split_objs[[i]])
  split_objs[[i]]@meta.data$celltype_coarse_AMfibro_merged <- original_meta[cells, , drop = FALSE]
  
  # 👇 ベクトルにしてから代入
  Idents(split_objs[[i]]) <- split_objs[[i]]@meta.data$celltype_coarse_AMfibro_merged[[1]]
}

# 各日付（day_group）で分割し connectome を再作成（新しい細胞分類で）
connectome_by_day_merged <- lapply(split_objs, function(seurat_obj) {
  CreateConnectome(
    seurat_obj,
    species = "mouse",
    p.values = TRUE,
    assay = "RNA",
    slot = "data",
    group.by = "celltype_coarse_AMfibro_merged"
  )
})



# 送信元: Mac_AM_group、受信側: Fibroblast_group
source_cells <- c("Mac_AM_group", "Mo_Ly6Chi", "DC_cDC2", "Epi_AT2", "Epi_AT1", 
                  "Endo_capillary", "Tcell_CD4T", "Tcell_CD8T", "Neutro", "NK")
target_cell <- "Fibroblast_group"

# 各日ごとに CircosPlot 描画
conflicted::conflicts_prefer(pheatmap::pheatmap)

for (day in names(connectome_by_day_merged)) {
  message("Plotting Circos: ", day)
  
  CircosPlot(connectome_by_day_merged[[day]],
             sources.include = source_cells,
             targets.include = target_cell,
             min.pct = 0.1,
             min.z = 1,
             lab.cex = 0.8,           # セクターラベル文字サイズ
             title = paste("Connectome Circos Plot -", day))
}


#LRP1だけに絞った解析
source_cells <- c("Mac_AM_group", "Mo_Ly6Chi", "DC_cDC2", "Epi_AT2", "Epi_AT1", 
                  "Endo_capillary", "Tcell_CD4T", "Tcell_CD8T", "Neutro", "NK")
target_cell <- "Fibroblast_group"

for (day in names(connectome_by_day_merged)) {
  message("Plotting Circos for day: ", day)
  
  conn <- connectome_by_day_merged[[day]]
  
  # 受容体が Lrp1、送信元が指定の細胞群、受信側がFibroblast_groupに限定
  conn_lrp1 <- conn %>%
    dplyr::filter(
      receptor == "Lrp1",
      source %in% source_cells,
      target == target_cell,
      percent.source > 0.1,
      percent.target > 0.1,
      weight_norm >= 1,
    )
  
  if (nrow(conn_lrp1) == 0) {
    message(" → Skipping ", day, " (no Lrp1 interactions in ", target_cell, ")")
    next
  }
  
  CircosPlot(conn_lrp1,
             lab.cex = 1.5,
             title = paste("LRP1-Fibroblast_group Connectome -", day))
}


#Figure3D
# condition から dayXX の "XX" だけ抽出
data_BLM$day_group2 <- gsub("day|BLM.*", "", data_BLM$condition)

# 数値として扱う
data_BLM$day_group2 <- as.numeric(data_BLM$day_group2)

# condition から最初の数字(1〜2桁)を抽出
data_BLM$day_group_num <- as.numeric(sub(".*day([0-9]+).*", "\\1", data_BLM$condition))

# 採用する7つの時系列だけ factor として使う
data_BLM$timepoint <- factor(
  data_BLM$day_group_num,
  levels = c(0, 3, 7, 14, 28, 42, 63),
  labels = c("day0", "day3", "day7", "day14", "day28", "day42", "day63")
)


fibro_only <- subset(data_BLM, subset = makoto.figure3 == "Fibroblast")

expr_df <- FetchData(fibro_only, vars = c("orig.ident", "timepoint", "Lrp1"))

# マウスごとの平均発現
mouse_summary <- expr_df %>%
  group_by(timepoint, orig.ident) %>%
  summarise(Lrp1_mean = mean(Lrp1, na.rm = TRUE)) %>%
  ungroup()

# timepoint ごとの平均（棒グラフ用）
time_summary <- mouse_summary %>%
  group_by(timepoint) %>%
  summarise(mean_Lrp1 = mean(Lrp1_mean))

library(ggplot2)

ggplot() +
  geom_col(data = time_summary,
           aes(x = timepoint, y = mean_Lrp1),
           fill = "gray70", width = 0.6) +
  geom_point(data = mouse_summary,
             aes(x = timepoint, y = Lrp1_mean),
             size = 3, alpha = 0.8) +
  theme_minimal() +
  ggtitle("LRP1 Expression in Fibroblasts (Per Mouse)") +
  ylab("Average LRP1 Expression") +
  xlab("Timepoint") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# マウスごとの平均
mouse_summary <- expr_df %>%
  group_by(orig.ident, timepoint) %>%
  summarise(Lrp1_mean = mean(Lrp1, na.rm = TRUE)) %>%
  ungroup()

# Excelに出力
write_xlsx(mouse_summary, "Fibroblast_LRP1_MouseLevel.xlsx")


#Figure3E

# Fibroblast だけを抽出（必要ならどちらか片方に絞ってください）
fibro <- subset(
  data_BLM,
  subset = makoto.figure3 %in% c("Fibroblast", "Fibroblast_Peribronchiolar")
)

# 時系列情報を Idents にセット
Idents(fibro) <- fibro$day_group
table(Idents(fibro))

# day42 Fibroblast で上昇している遺伝子を抽出
deg_day42 <- FindMarkers(
  fibro,
  ident.1 = "day42",     # 比較したい群
  ident.2 = NULL,        # 他の全 timepoint と比較
  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE        # day42 で高い遺伝子だけ
)

# 行名を列に
deg_day42 <- deg_day42 %>%
  mutate(gene = rownames(deg_day42))

enriched_day42 <- deg_day42 %>%
  filter(p_val_adj < 0.05) %>%      # 多重検定後有意
  arrange(desc(avg_log2FC))        # log2FC の高い順

# 必要ならトップ100だけ
# enriched_day42 <- enriched_day42 %>% slice_max(avg_log2FC, n = 100)
write_xlsx(enriched_day42,
           path = "Fibroblast_day42_enriched_genes.xlsx")






