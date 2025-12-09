#IP303例データ解析

setwd("C:\\Users\\myama\\Desktop\\IP303例")
data_IP<- readRDS("C:\\Users\\myama\\Desktop\\IP303例\\HC_NSIP_ILD_HP_Myo_SSc_303samples.rds")


# biomaRtを使用してENSEMBL IDを遺伝子シンボルに変換
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 遺伝子名（ENSEMBL ID）を取得
ensembl_ids <- rownames(GetAssayData(data_IP, slot = "counts"))

# ENSEMBL IDを遺伝子シンボルに変換
gene_conversion <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# ENSEMBL IDを遺伝子シンボルにマッピング
gene_conversion <- gene_conversion[gene_conversion$hgnc_symbol != "", ]
gene_mapping <- setNames(gene_conversion$hgnc_symbol, gene_conversion$ensembl_gene_id)

# Seuratオブジェクト内の遺伝子名を遺伝子シンボルに更新
counts_data <- GetAssayData(data_IP, slot = "counts")
new_rownames <- ifelse(!is.na(gene_mapping[rownames(counts_data)]), gene_mapping[rownames(counts_data)], rownames(counts_data))
unique_new_rownames <- make.unique(new_rownames)
rownames(counts_data) <- unique_new_rownames

# 'data'スロットも同様に更新
data_data <- GetAssayData(data_IP, slot = "data")
rownames(data_data) <- unique_new_rownames

# 元のオブジェクトからAssayを作成し、それぞれ設定
assay <- CreateAssayObject(counts = counts_data)
assay <- SetAssayData(assay, slot = "data", new.data = data_data)

# Seuratオブジェクトに設定
data_IP[['RNA']] <- assay

# 確認
head(rownames(GetAssayData(data_IP, slot = "counts")))


DimPlot(data_IP, reduction = "umap")
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_1", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_2", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_4", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_5", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "original_ann_level_1", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "original_ann_level_2", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "original_ann_level_3", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "original_ann_level_4", label =TRUE)
DimPlot(data_IP, reduction = "umap", group.by = "original_ann_level_5", label =TRUE)

DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = TRUE, repel = TRUE) + theme_void() + NoLegend() + theme(plot.margin = margin(10, 10, 10, 10), text = element_text(size = 18))
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = TRUE, repel = TRUE) + scale_color_viridis_d(option = "mako") + theme_void() + NoLegend() + theme(plot.margin = margin(5, 5, 5, 5), text = element_text(size = 16))
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = TRUE, repel = TRUE) + scale_color_viridis_d(option = "cividis") + theme_void() + NoLegend() + theme(plot.margin = margin(5, 5, 5, 5), text = element_text(size = 16))
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = TRUE, repel = TRUE) + scale_color_viridis_d(option = "inferno") + theme_void() + NoLegend() + theme(plot.margin = margin(5, 5, 5, 5), text = element_text(size = 16))
DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = TRUE, repel = TRUE) + scale_color_viridis_d(option = "turbo") + theme_void() + NoLegend() + theme(plot.margin = margin(5, 5, 5, 5), text = element_text(size = 16))


Idents(data_IP) <- "ann_level_3"
n_clust <- length(levels(data_IP))

# 落ち着いたブルー系パレットを作成
blue_pal <- colorRampPalette(brewer.pal(9, "Blues"))(n_clust)

p_blue <- DimPlot(
  data_IP,
  reduction = "umap",
  label = TRUE,           # ラベルは今まで通り
  repel = TRUE,
  cols  = blue_pal,
  raster = TRUE           # 点が多いので raster=TRUE のままでOK
) +
  theme_void() +          # 軸線・目盛・文字を全部消す
  theme(
    plot.title = element_text(hjust = 0, size = 18),
    legend.position = "none"
  )

p_blue


Idents(data_IP) <- "ann_level_3"
n_clust <- length(levels(data_IP))

# ブラウン系（少しグリーンも入るが落ち着いたトーン）
brown_pal <- colorRampPalette(brewer.pal(11, "BrBG"))(n_clust)

p_brown <- DimPlot(
  data_IP,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  cols  = brown_pal,
  raster = TRUE
) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 18),
    legend.position = "none"
  )

p_brown

#Figure用分類作成
# ベースは ann_level_3
data_IP$makoto.figure1 <- as.character(data_IP$ann_level_3)

ann4 <- data_IP$original_ann_level_4

## --- Alveolar macrophages ---
alveolar_labels <- c(
  "Alveolar macrophages",
  "Alveolar macrophages MARCO-",
  "Alveolar macrophages MARCO+",   # あれば
  "Alveolar macrophages SPP1+"     # あれば
)

data_IP$makoto.figure1[ann4 %in% alveolar_labels] <- "Alveolar macrophages"

## --- Interstitial macrophages ---
interstitial_labels <- c(
  "Interstitial macrophages",
  "Macrophages SPP1 high",
  "Mac proliferating",
  "Macrophages MARCO-",
  "Macrophages MARCO+"             # あれば
)

data_IP$makoto.figure1[ann4 %in% interstitial_labels] <- "Interstitial macrophages"

## --- 血中モノサイト 3 サブセット ---
data_IP$makoto.figure1[ann4 %in% "Classical monocytes"]   <- "Classical monocytes"
data_IP$makoto.figure1[ann4 %in% "Intermediate monocytes"]<- "Intermediate monocytes"
data_IP$makoto.figure1[ann4 %in% "Non-classical monocytes"] <- "Non-classical monocytes"

## DC monocyte など半端なものは必要に応じてどこかへ
data_IP$makoto.figure1[ann4 %in% "DC monocyte"] <- "Classical monocytes"  # 例

# 好きな順番で level を整理
macro_levels <- c(
  "Alveolar macrophages",
  "Interstitial macrophages",
  "Classical monocytes",
  "Intermediate monocytes",
  "Non-classical monocytes"
)

conflicts_prefer(GSEABase::setdiff)
conflicts_prefer(lubridate::setdiff)
conflicts_prefer(GenomicRanges::setdiff)
conflicts_prefer(base::setdiff)

other_levels <- setdiff(unique(data_IP$makoto.figure1), macro_levels)

data_IP$makoto.figure1 <- factor(
  data_IP$makoto.figure1,
  levels = c(macro_levels, sort(other_levels))
)

# 確認
table(data_IP$makoto.figure1)

Idents(data_IP) <- "makoto.figure1"
n_clust <- length(levels(data_IP))

blue_pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n_clust)

p_fig1 <- DimPlot(
  data_IP,
  reduction = "umap",
  label = TRUE,
  repel  = TRUE,
  cols   = blue_pal,
  raster = TRUE
) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 18),
    legend.position = "none"
  ) +
  ggtitle("makoto.figure1")

p_fig1

p_fig1_no_raster <- p_fig1 + theme_void() + guides(colour = "none")
tiff("figure_no_raster.tiff", width=10, height=10, units="in", res=1200, type="cairo")
print(p_fig1_no_raster)
dev.off()

#グレースケール
p_lrp1_black <- FeaturePlot(
  data_IP,
  features   = "LRP1",
  min.cutoff = "q60",         # 低発現はほぼ真っ白に飛ばす
  max.cutoff = "q99",         # ごく高発現だけ黒に
  pt.size    = 0.2,
  order      = TRUE,          # 高発現ポイントを上に描画
  raster     = FALSE          # ベクター出力（高画質）
) +
  # cols は使わず、scale_color_gradientn でまとめて指定
  scale_color_gradientn(
    colors = c("grey95", "grey70", "black")  # コントラスト強め
  ) +
  theme_void() +               # 軸線・数字・枠など完全削除
  labs(title = "Lrp1 Expression in Fibroblasts") +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 16),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

# 画面で確認
p_lrp1_black

# 凡例なしバージョン
p_lrp1_noleg <- p_lrp1_black + guides(colour = "none")

# ---- 高画質 TIFF 出力 ----
tiff("LRP1_fibro_strongContrast.tiff",
     width = 10, height = 10, units = "in", res = 1200, type = "cairo")
print(p_lrp1_noleg)
dev.off()



# UMAPプロットを作成し、カスタムテーマを適用
umap_plot <- DimPlot(data_IP, reduction = "umap", group.by = "ann_level_3", label = FALSE, pt.size = 1.0) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),  # 枠線を消す
    panel.grid.major = element_blank(),  # 主なグリッド線を消す
    panel.grid.minor = element_blank(),  # 小さなグリッド線を消す
    panel.background = element_blank(),  # 背景を消す
    axis.line = element_blank(),  # 軸の線を消す
    axis.ticks = element_blank(),  # 軸の目盛りを消す
    axis.text = element_blank(),  # 軸のテキストを消す
    axis.title = element_blank(),  # 軸のタイトルを消す
    legend.position = "none"  # 凡例を消す
  )

print(umap_plot)


# `makoto.figure1`がNAでないデータのみを使用
valid_cells <- rownames(data_IP@meta.data)[!is.na(data_IP@meta.data$makoto.figure1)]
data_subset <- subset(data_IP, cells = valid_cells)

# 必要な細胞種のみを選択するリスト
cell_order <- c("B cell lineage", "T cell lineage", "Dendritic cells", "EC capillary", "AT1", "AT2", 
                "Fibroblasts", "Myofibroblasts", 
                "Classical monocytes", "Intermediate monocytes", "Non-classical monocytes", "Alveolar macrophages", "Interstitial macrophages", "Vascular smooth muscle")

# 指定した細胞種のみを含むデータを抽出
data_filtered <- subset(data_subset, subset = makoto.figure1 %in% cell_order)


# 行名の重複がないか確認し、修正する
rownames(data_filtered@meta.data) <- make.unique(as.character(rownames(data_filtered@meta.data)))
rownames(data_filtered) <- make.unique(as.character(rownames(data_filtered)))

# Dot plotの作成
dot_plot <- DotPlot(
  object = data_filtered,
  features = "LRP1",
  group.by = "makoto.figure1",  # クラスターでグループ化
  cols = c("lightgrey", "black"),  # カラースケールの設定
  dot.scale = 6  # ドットサイズのスケール
) + 
  labs(
    title = "LRP1 Expression Across makoto.figure1",
    x = "Cell Types",
    y = "Expression Level"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # X軸のラベルを45度回転
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  # 格子を消す
  )

# プロットの表示
print(dot_plot)

dot_plot2 <- dot_plot +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = "LRP1_dotplot_white_bg_1200dpi.tiff",
  plot     = dot_plot2,
  device   = "tiff",
  dpi      = 1200,
  width    = 10,
  height   = 8,
  units    = "in",
  bg       = "white",      # ← 背景を確実に白にする
  compression = "lzw"
)


# LRP1発現をcelltype毎にバイオリンプロットで表示、格子を消す
# Seuratオブジェクトの条件グループを設定 Healthyとそれ以外のILDの2郡
data_IP$condition_group <- ifelse(data_IP$lung_condition == "Healthy", "Healthy", "ILD")
# 必要な細胞種のみを選択するリスト
cell_order <- c("B cell lineage", "T cell lineage", "Dendritic cells", "EC capillary", "AT1", "AT2", 
                "Fibroblasts", "Myofibroblasts", 
                "Classical monocytes", "Intermediate monocytes", "Non-classical monocytes", "Alveolar macrophages", "Interstitial macrophages", "Vascular smooth muscle")

# 必要な細胞種のみを含むデータを抽出
data_filtered <- subset(data_IP, subset = makoto.figure1 %in% cell_order)
table(data_IP$lung_condition, data_IP$condition_group)

conflicts_prefer(matrixStats::count)
conflicts_prefer(dplyr::count)
df <- FetchData(data_filtered, vars = c("LRP1", "makoto.figure1", "condition_group"))
counts_table <- df %>%
  count(makoto.figure1, condition_group) %>%
  pivot_wider(names_from = condition_group,
              values_from = n,
              values_fill = 0)

counts_table

# バイオリンプロットの作成
vln_plot <- VlnPlot(
  object = data_filtered, 
  features = "LRP1", 
  group.by = "makoto.figure1", 
  split.by = "condition_group", 
  split.plot = TRUE, 
  pt.size = 0  # ドットを非表示にする
) +
  ggtitle("LRP1 Expression Across Cell Types (Healthy vs ILD)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  # 格子を消す
  ) +
  scale_fill_manual(values = c("Healthy" = "#FCFCFC", "ILD" = "#828282"))

# プロットを表示
print(vln_plot)

vln_plot2 <- vln_plot +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = "LRP1_vlnplot_white_bg_1200dpi.tiff",
  plot     = vln_plot2,
  device   = "tiff",
  dpi      = 1200,
  width    = 10,
  height   = 8,
  units    = "in",
  bg       = "white",      # ← 背景を確実に白にする
  compression = "lzw"
)

# 必要なメタデータと LRP1 発現を取り出す
df <- FetchData(data_filtered, vars = c("LRP1", "makoto.figure1", "condition_group"))

results <- df %>%
  group_by(makoto.figure1) %>%
  summarise(
    n_Healthy = sum(condition_group == "Healthy"),
    n_ILD     = sum(condition_group == "ILD"),
    p_value   = ifelse(
      n_Healthy > 0 & n_ILD > 0,                 # 両群とも少なくとも1細胞あるときだけ
      wilcox.test(LRP1 ~ condition_group)$p.value,
      NA_real_                                   # どちらか0なら p値は NA
    ),
    .groups = "drop"
  )

# FDR補正（BH）
results$FDR <- p.adjust(results$p_value, method = "BH")
results


# FibroblastsとMyofibroblastsのみを抽出
Idents(data_IP) <- data_IP@meta.data$original_ann_level_3
fibroblasts_myofibroblasts_data <- subset(data_IP, idents = c("Fibroblasts", "Myofibroblasts"))

DimPlot(fibroblasts_myofibroblasts_data, reduction = "umap", group.by = "ann_level_4", label =TRUE)
DimPlot(fibroblasts_myofibroblasts_data, reduction = "umap", group.by = "original_ann_level_4", label =TRUE)

# SCTransformによるデータの正規化と再クラスタリング
plan("sequential")
options(future.globals.maxSize = 20 * 1024^3)  # 20GB くらい
fibroblasts_myofibroblasts_data <- SCTransform(fibroblasts_myofibroblasts_data, verbose = TRUE)
fibroblasts_myofibroblasts_data <- RunPCA(fibroblasts_myofibroblasts_data, verbose = TRUE)
fibroblasts_myofibroblasts_data <- FindNeighbors(fibroblasts_myofibroblasts_data, dims = 1:30, verbose = TRUE)
fibroblasts_myofibroblasts_data <- FindClusters(fibroblasts_myofibroblasts_data, resolution = 0.5, verbose = TRUE)
fibroblasts_myofibroblasts_data <- RunUMAP(fibroblasts_myofibroblasts_data, dims = 1:30, verbose = TRUE)
DimPlot(fibroblasts_myofibroblasts_data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(fibroblasts_myofibroblasts_data, features = c("LRP1"), max.cutoff = 3, cols = c("#FFE4C4", "darkred", pt.size = 1.5))


# FeaturePlotを作成し、カスタムテーマを適用
feature_plot <- FeaturePlot(fibroblasts_myofibroblasts_data, features = c("LRP1"), max.cutoff = 3, 
                            cols = c("#FFE4C4", "darkred"), pt.size = 1.5) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),  # 枠線を消す
    panel.grid.major = element_blank(),  # 主なグリッド線を消す
    panel.grid.minor = element_blank(),  # 小さなグリッド線を消す
    panel.background = element_blank(),  # 背景を消す
    axis.line = element_blank(),  # 軸の線を消す
    axis.ticks = element_blank(),  # 軸の目盛りを消す
    axis.text = element_blank(),  # 軸のテキストを消す
    axis.title = element_blank()  # 軸のタイトルを消す
  ) +
  
  scale_color_gradientn(
    colors = c("#FFE4C4", "darkred"),  # カラースケールの始点と終点の色
    limits = c(0, 3),  # カラースケールのデータ範囲
    values = c(0, 0.1, 0.5, 1)  # 色の割り当てを制御するポイント
  )

print(feature_plot)

#fibroblastsにおけるクラスター毎のLRP1発現の箱ひげ図
#SCTよりRNAアッセイのほうが滑らかなdotになる
data <- fibroblasts_myofibroblasts_data@meta.data %>%
  mutate(LRP1_expression = fibroblasts_myofibroblasts_data@assays$SCT@data["LRP1", ]) %>%
  filter(!is.na(LRP1_expression))  # 発現がNAでないデータのみを使用

ggplot(data, aes(x = factor(seurat_clusters), y = LRP1_expression)) +
  geom_boxplot(outlier.size = 0, fill = "grey", color = "black") +  # 箱ひげ図
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkred") +        # ドットプロット
  labs(title = "LRP1 Expression Across Clusters", 
       x = "Clusters", 
       y = "LRP1 Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # X軸ラベルの調整
    panel.grid.major = element_blank(),  # 主なグリッド線を消す
    panel.grid.minor = element_blank()   # 小さなグリッド線を消す
  )

#fibroblastsにおけるクラスター毎のLRP1発現の有意差検定
data <- fibroblasts_myofibroblasts_data@meta.data %>%
  mutate(Lrp1_expression = fibroblasts_myofibroblasts_data@assays$SCT@data["LRP1", ]) %>%
  dplyr::filter(!is.na(Lrp1_expression))
kruskal_test <- kruskal.test(Lrp1_expression ~ seurat_clusters, data = data)
dunn_test <- dunnTest(Lrp1_expression ~ seurat_clusters, data = data, method = "bonferroni")
dunn_results <- as.data.frame(dunn_test$res)
significant_pairs <- dunn_results %>%
  dplyr::filter(P.adj < 0.05) %>%
  dplyr::select(Comparison, P.adj)
print(dunn_results)
print(significant_pairs)
cluster_11_results <- dunn_results %>%
  dplyr::filter(grepl("11 -", Comparison) | grepl(" - 11", Comparison))
print(cluster_11_results)


#線維芽細胞のtsukui論文参考のアノテーション


# === 線維芽細胞サブセットの抽出 ========================================

fibro_targets <- c("Adventitial fibroblasts", "Alveolar fibroblasts", 
                   "Fibroblasts PLIN2+", "Fibroblasts WIF1+ CHRM2+", 
                   "Peribronchial fibroblasts", "Nerve-associated fibroblasts",
                   "Subpleural fibroblasts", "Secondary crest myofibroblast",
                   "Activated myofibroblasts", "Lipofibroblasts")

fibro_subset <- subset(data_IP, subset = original_ann_level_4 %in% fibro_targets)


# === 通常の前処理（Normalize → PCA → UMAP → クラスタリング） =======================
fibro_subset <- NormalizeData(fibro_subset)
fibro_subset <- FindVariableFeatures(fibro_subset)
fibro_subset <- ScaleData(fibro_subset)
fibro_subset <- RunPCA(fibro_subset)
fibro_subset <- RunUMAP(fibro_subset, dims = 1:15)
fibro_subset <- FindNeighbors(fibro_subset, dims = 1:15)
fibro_subset <- FindClusters(fibro_subset, resolution = 0.5)

marker_genes <- c("TCF21", "PDGFRA", "SAA3", "LCN2", "CCL2", "CXCL14", 
                  "COL1A1", "CTHRC1", "POSTN", "HSPA1A", "HSPA1B", 
                  "MKI67", "TOP2A", "HHIP", "PI16", "DCN", "NPNT", "CES1D", "WNT2A", "CDKN1A", "GDF15", "WNT5A", "CDC20", "CDCA8", "SFRP4", "IFI27", "CXCL12", "CXCL2", "SFRP2", "INMT")

DotPlot(fibro_subset, features = marker_genes) + RotatedAxis()

#上記で分類した０～１６のSuerat Clustersを再度data_IPにメタデータとして記録する
data_IP$original_ann_level_4_seurat_clusters <- as.character(data_IP$original_ann_level_4)
fibro_cluster_ids <- as.character(Idents(fibro_subset))
names(fibro_cluster_ids) <- colnames(fibro_subset)
data_IP$original_ann_level_4_seurat_clusters[names(fibro_cluster_ids)] <- fibro_cluster_ids

# 再アノテーション：クラスタに任意の名前をつける
refined_labels <- c(
  "0" = "Adventitial fibroblasts",
  "1" = "Alveolar fibroblasts",
  "2" = "Inflammatory fibroblasts2",
  "3" = "Adventitial fibroblasts",  # 評価困難なため保留（除外推奨）
  "4" = "Fibrotic fibroblasts",
  "5" = "Fibrotic fibroblasts",
  "6" = "Fibrotic fibroblasts",
  "7" = "Inflammatory fibroblasts1",
  "8" = "Adventitial fibroblasts",
  "9" = "Inflammatory fibroblasts1",
  "10" = "Fibrotic fibroblasts",
  "11" = "Adventitial fibroblasts",
  "12" = "Adventitial fibroblasts",
  "13" = "Peribronchial fibroblasts",
  "14" = "Alveolar fibroblasts",
  "15" = "Alveolar fibroblasts",
  "16" = "Fibrotic fibroblasts"
)

# クラスタIDに対応する注釈ベクトルを作成
refined_labels <- refined_labels[as.character(Idents(fibro_subset))]
names(refined_labels) <- colnames(fibro_subset)
fibro_subset$refined_fibro_anno <- refined_labels

# === data_IP に refined 注釈を統合する =====================================

# 1. 新しいメタデータ列として初期化（元の注釈をベースに）
data_IP$original_ann_level_4_refined <- as.character(data_IP$original_ann_level_4)

# 2. fibro_subset の refined 注釈を抽出
fibro_labels <- fibro_subset$refined_fibro_anno
names(fibro_labels) <- colnames(fibro_subset)

# 3. 対象細胞のアノテーションを上書き
data_IP$original_ann_level_4_refined[names(fibro_labels)] <- fibro_labels

# === 確認 ================================================
# 該当細胞の再注釈されたカテゴリ数を確認
table(data_IP$original_ann_level_4_refined[names(fibro_labels)])

# 全体にどのラベルが含まれているか確認
unique(data_IP$original_ann_level_4_refined)

# NAを除くセルの名前だけを取得
valid_cells <- colnames(data_IP)[!is.na(data_IP$original_ann_level_4_refined)]

# それらの細胞だけを取り出す
data_IP_clean <- subset(data_IP, cells = valid_cells)
DotPlot(data_IP_clean, features = "LRP1", group.by = "original_ann_level_4_refined") + RotatedAxis()

#=====ここまでで　つくい論文の線維芽細胞アノテーション完了、original_ann_level_refinedのmeta.data作成した====

# 表示したい細胞型（順番：Fibroblast → Macrophage）
ordered_celltypes <- c(
  # Fibroblast系
  "Adventitial fibroblasts",
  "Alveolar fibroblasts",
  "Inflammatory fibroblasts1",
  "Inflammatory fibroblasts2",
  "Fibrotic fibroblasts",
  "Peribronchial fibroblasts",
  # Macrophage系
  "Alveolar macrophages",
  "Interstitial macrophages",
  "Macrophages MARCO-",
  "Macrophages SPP1 high"
)

# セルを抽出
valid_cells <- WhichCells(data_IP, expression = original_ann_level_4_refined %in% ordered_celltypes)
data_IP_clean <- subset(data_IP, cells = valid_cells)

# 図示前に factor の順番を設定（並び順調整）
data_IP_clean$original_ann_level_4_refined <- factor(
  data_IP_clean$original_ann_level_4_refined,
  levels = ordered_celltypes
)

# DotPlot 作成（グレースケール）
DotPlot(data_IP_clean,
        features = "LRP1",
        group.by = "original_ann_level_4_refined",
        cols = c("grey90", "black")) +
  RotatedAxis() +
  ggtitle("LRP1 expression in fibroblast and macrophage subtypes") +
  theme(plot.title = element_text(hjust = 0.5))

#=====つくい論文の線維芽細胞4種の時系列解析=================
#Monocle３解析

target_fibrotypes <- c("Alveolar fibroblasts", "Inflammatory fibroblasts1", 
                       "Inflammatory fibroblasts2", "Fibrotic fibroblasts")

# サブセット抽出
fib_subset <- subset(data_IP, cells = WhichCells(data_IP, expression = original_ann_level_4_refined %in% target_fibrotypes))

# 発現行列、細胞情報、遺伝子情報を抽出
expr_matrix <- GetAssayData(fib_subset, slot = "counts")
cell_metadata <- fib_subset@meta.data
gene_metadata <- data.frame(
  gene_short_name = rownames(expr_matrix),
  row.names = rownames(expr_matrix)
)

# Monocle3オブジェクトを作成
fibroblast_cds <- new_cell_data_set(
  expression_data = expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

colData(fibroblast_cds)$cell_type <- fib_subset$original_ann_level_4_refined
fibroblast_cds <- preprocess_cds(fibroblast_cds, num_dim = 50)
fibroblast_cds <- reduce_dimension(fibroblast_cds, reduction_method = "UMAP")
fibroblast_cds <- cluster_cells(fibroblast_cds, reduction_method = "UMAP")
fibroblast_cds <- learn_graph(fibroblast_cds)
plot_cells(fibroblast_cds, color_cells_by = "cell_type")

# ルート状態を設定（: Alveolar　Fibroblasts）
root_cells <- colnames(fibroblast_cds)[colData(fibroblast_cds)$cell_type == "Alveolar fibroblasts"]
fibroblast_cds <- order_cells(fibroblast_cds, root_cells = root_cells)

# 擬似時系列プロット
plot_cells(
  fibroblast_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = FALSE, # 遷移グラフを非表示にする場合
  cell_size = 1.0 # ドットサイズを指定（デフォルトは1）
) +
  theme_void() + # 軸、ラベル、数字をすべて非表示に
  theme(legend.title = element_text(size = 14), # 凡例タイトルのサイズを調整（任意）
        legend.text = element_text(size = 12)) # 凡例ラベルのサイズを調整（任意）

# LRP1の発現可視化
plot_cells(
  fibroblast_cds,
  genes = "LRP1",
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 1.0  # ドットを大きく
) +
  theme_void() +  # 軸や目盛りを非表示に
  ggtitle("LRP1 expression in trajectory") +
  theme(
    plot.title = element_text(size = 16, face = "bold")  # タイトルの調整（任意）
  )

# 細胞種ごとの可視化（original_ann_level_4_refined）
plot_cells(
  fibroblast_cds,
  color_cells_by = "original_ann_level_4_refined",
  show_trajectory_graph = FALSE, # 遷移グラフを非表示にする場合
  cell_size = 1.0 # ドットサイズを指定（デフォルトは1）
) +
  theme_void() + # 軸、ラベル、数字をすべて非表示に
  theme(legend.title = element_text(size = 14), # 凡例タイトルのサイズを調整（任意）
        legend.text = element_text(size = 12)) # 凡例ラベルのサイズを調整（任意）

#時系列で遺伝子発現の推移を見るグラフの作図
# 特定の遺伝子の発現データを抽出
genes_of_interest <- c("LRP1", "INMT", "TCF21", "CXCL12", "CXCL14", "CTHRC1", "COL1A1")

# Monocle3でpseudotimeに沿った発現をプロット
monocle3::plot_genes_in_pseudotime(
  fibroblast_cds[genes_of_interest, ]
)

#monocle3のデータをSCトランスフォーム後のUMAPで表示させるスクリプト
options(future.globals.maxSize = 2 * 1024^3)  # 2GB
fib_subset <- SCTransform(fib_subset, verbose = FALSE)
fib_subset <- SCTransform(fib_subset, verbose = FALSE)
fib_subset <- RunPCA(fib_subset)
fib_subset <- RunUMAP(fib_subset, dims = 1:20)
# pseudotime を計算済みの Monocle3 オブジェクト
pseudotime_vec <- pseudotime(fibroblast_cds)
umap_coords <- Embeddings(fib_subset, reduction = "umap")
#UMAP座標とPsuedtimeを統合
plot_df <- data.frame(
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Pseudotime = pseudotime_vec[rownames(umap_coords)]
)

ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Pseudotime)) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "magma") +
  theme_void() +
  labs(title = "Pseudotime projected onto Seurat UMAP") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# 細胞種をIdentsに指定
Idents(fib_subset) <- fib_subset$original_ann_level_4_refined

# UMAP上に細胞種を表示
DimPlot(
  fib_subset,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  pt.size = 1.2 # ドットサイズを大きく
) +
  ggtitle("Cell types on UMAP") +
  theme_void() +  # 軸、目盛り、枠線を非表示
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

FeaturePlot(
  fib_subset,
  features = "LRP1",
  reduction = "umap",
  pt.size = 1.0
) +
  ggtitle("LRP1 expression on UMAP") +
  theme_void() +  # 軸、目盛り、枠線を非表示
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

FeaturePlot(
  fib_subset,
  features = "LRP1",
  reduction = "umap",
  pt.size = 1.0
) +
  ggtitle("LRP1 expression on UMAP") +
  theme_void() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_color_gradientn(
    colours = c("grey90", "orange", "red"),
    values = scales::rescale(c(0, 2, 3)),  # 値1〜3付近で急激に色変化
    oob = scales::squish
  )


#CytoTRACEによる時系列解析
# CytoTRACE解析のための準備
expression_matrix <- as.matrix(fib_subset@assays$RNA@counts)

# CytoTRACEの実行
cytotrace_results <- CytoTRACE(expression_matrix)
cytotrace_scores <- cytotrace_results$CytoTRACE
cell_names <- colnames(fib_subset)
names(cytotrace_scores) <- cell_names
fib_subset$CytoTRACE <- cytotrace_scores
colData(fibroblast_cds)$CytoTRACE <- cytotrace_scores

# CytoTRACEスコアを色分けに用いたPseudotimeプロットの描画
plot_cells(
  fibroblast_cds,
  color_cells_by = "CytoTRACE",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_branch_points = TRUE,
  cell_size = 1.0
) +
  scale_color_gradient(low = "#FFE4C4", high = "darkgreen") +
  theme_minimal() +
  labs(title = "CytoTRACE Scores on Trajectory") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


FeaturePlot(fib_subset, features = "CytoTRACE", reduction = "umap", cols = c("lightgrey", "darkgreen"), pt.size = 1.2) +
  scale_color_gradient(low = "lightgrey", high = "darkgreen") +
  theme_minimal() +
  labs(title = "CytoTRACE Scores on Pseudotime Trajectory") +
  theme(
    panel.grid.major = element_blank(),  # メジャーグリッドラインを非表示
    panel.grid.minor = element_blank()   # マイナーグリッドラインを非表示
  )

DimPlot(
  fib_subset,
  group.by = "original_ann_level_4_seurat_clusters",
  reduction = "umap",
  label = TRUE,                 # 各クラスタにラベルを表示
  label.size = 4,
  repel = TRUE,                 # ラベルが重ならないように調整
  cols = NULL                  # 自動で色分け（必要ならカスタム色も可）
) +
  theme_minimal() +
  labs(title = "UMAP by original_ann_level_4_seurat_clusters") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


# クラスタ番号は original_ann_level_4_seurat_clusters に記録されている前提
cluster_cytotrace_summary <- fib_subset@meta.data %>%
  group_by(original_ann_level_4_seurat_clusters) %>%
  summarise(mean_CytoTRACE = mean(CytoTRACE, na.rm = TRUE)) %>%
  arrange(desc(mean_CytoTRACE))  # 高いほど未分化＝始点側

print(cluster_cytotrace_summary)

# pseudotime値をメタデータに追加
fibroblast_cds$pseudotime <- monocle3::pseudotime(fibroblast_cds)

# クラスタごとの平均 pseudotime を計算
pseudotime_summary <- as.data.frame(colData(fibroblast_cds)) %>%
  group_by(original_ann_level_4_seurat_clusters) %>%
  summarise(mean_pseudotime = mean(pseudotime, na.rm = TRUE)) %>%
  arrange(mean_pseudotime)

print(pseudotime_summary)

# 遺伝子リスト
genes_of_interest <- c("LRP1", "INMT", "TCF21", "CXCL12", "CXCL14", "ACTA2", "COL1A1", "COL3A1", "FN1", "COL1A2")


# メタデータと発現を結合（fib_subsetはSeuratオブジェクト）
plot_data <- FetchData(fib_subset, vars = c(genes_of_interest, "CytoTRACE")) %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "gene", values_to = "expression")

# ログスケール用（0を避けるため+1e-3など追加）
plot_data$expression_log <- log10(plot_data$expression + 1e-3)

# geneごとのプロットを作成
gene_plots <- lapply(genes_of_interest, function(g) {
  plot_data %>%
    dplyr::filter(gene == g) %>%
    ggplot(aes(x = CytoTRACE, y = expression_log, color = CytoTRACE)) +
    geom_point(size = 0.4, alpha = 0.6) +
    geom_smooth(color = "black", method = "loess", se = FALSE) +
    scale_color_viridis_c() +
    labs(title = g, x = "CytoTRACE Score", y = "Log10(Expression)") +
    theme_minimal(base_size = 12)
})


# patchworkで並列配置
wrap_plots(gene_plots, ncol = 2)


# CytoTRACEをUMAP上にマッピング
umap_df <- as.data.frame(Embeddings(fib_subset, "umap"))
umap_df$CytoTRACE <- fib_subset$CytoTRACE

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_viridis_c(option = "D", limits = c(0, 1)) +
  labs(title = "CytoTRACE Score on UMAP", color = "CytoTRACE") +
  theme_minimal()


# 色スケールを同じにした version
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE)) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_viridis_c(option = "D", limits = c(0, 1)) +
  labs(title = "CytoTRACE Score on UMAP", color = "CytoTRACE") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),        # 格子線を消す
    axis.text = element_blank(),         # 軸の数字を消す
    axis.ticks = element_blank(),        # 軸目盛りを消す
    axis.title = element_blank(),        # 軸ラベル (UMAP_1, UMAP_2) を消す
    axis.line = element_blank()          # 軸線そのものを消す
  )




# クラスタと細胞種の対応を定義
cluster_to_celltype <- c(
  "1" = "Fibrotic fibroblasts",
  "5" = "Fibrotic fibroblasts",
  "9" = "Fibrotic fibroblasts",
  "16" = "Fibrotic fibroblasts",
  "7" = "Inflammatory fibroblasts2",
  "6" = "Inflammatory fibroblasts1",
  "2" = "Inflammatory fibroblasts1",
  "10" = "Alveolar fibroblasts",
  "4" = "Alveolar fibroblasts",
  "15" = "Alveolar fibroblasts",
  "14" = "Alveolar fibroblasts"
)

# 新しい列にアノテーションを追加（元のアノテーションが残るように）
fib_subset@meta.data$reannotated_fibro_type <- cluster_to_celltype[
  as.character(fib_subset@meta.data$original_ann_level_4_seurat_clusters)
]

# 確認
table(fib_subset@meta.data$reannotated_fibro_type, useNA = "always")

# 再アノテーション情報をベクトルとして取り出す
# 各オブジェクトの細胞名を取得
conflicted::conflicts_prefer(dplyr::desc)

# 各オブジェクトの細胞バーコードを取得
fib_barcodes <- colnames(fib_subset)
data_barcodes <- colnames(data_IP)

# クリーンなバーコードに変換（"_”以降が真のバーコードと仮定）
fib_clean_barcodes <- sub(".*_", "", fib_barcodes)
data_clean_barcodes <- sub(".*_", "", data_barcodes)

# 各オブジェクトにクリーンバーコードを追加
fib_subset@meta.data$clean_barcode <- fib_clean_barcodes
data_IP@meta.data$clean_barcode <- data_clean_barcodes

# fib_subset のアノテーション情報をデータフレーム化
anno_df <- data.frame(
  clean_barcode = fib_subset@meta.data$clean_barcode,
  new_annotation = fib_subset@meta.data$reannotated_fibro_type,
  row.names = NULL
)

# 重複除去
anno_df_unique <- anno_df[!duplicated(anno_df$clean_barcode), ]

# data_IPのmeta.dataの最低限の列を取り出し
meta_df_minimal <- data_IP@meta.data[, c("clean_barcode", "original_ann_level_4")]

# マージ（left joinでnew_annotationを追加）
merged_meta <- meta_df_minimal %>%
  left_join(anno_df_unique, by = "clean_barcode")

# 新しい注釈列 reannotated_fibro_type を作成
merged_meta$reannotated_fibro_type <- ifelse(
  is.na(merged_meta$new_annotation),
  as.character(merged_meta$original_ann_level_4),
  as.character(merged_meta$new_annotation)
)

# meta.dataに再格納
data_IP@meta.data$reannotated_fibro_type <- merged_meta$reannotated_fibro_type
rownames(data_IP@meta.data) <- rownames(data_IP@meta.data)  # 念のため元の順序保持

# ✔ 確認
table(data_IP@meta.data$reannotated_fibro_type, useNA = "always")


# UMAP座標をデータフレームに変換
umap_df <- as.data.frame(Embeddings(fib_subset, "umap"))
umap_df$reannotated_fibro_type <- fib_subset@meta.data$reannotated_fibro_type

# UMAPプロットを作成
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = reannotated_fibro_type)) +
  geom_point(size = 1.2, alpha = 0.8) +
  labs(title = "UMAP with Reannotated Fibroblast Types", color = "Cell Type") +
  scale_color_manual(values = c(
    "Fibrotic fibroblasts" = "firebrick",
    "Inflammatory fibroblasts1" = "steelblue",
    "Inflammatory fibroblasts2" = "skyblue",
    "Alveolar fibroblasts" = "forestgreen"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),         # 軸タイトル（UMAP_1, UMAP_2）削除
    axis.text = element_blank(),          # 軸目盛りのラベル削除
    axis.ticks = element_blank(),         # 軸目盛りの線削除
    axis.line = element_blank()           # 軸線削除（theme_minimalには元々ないが保険で）
  )

print(p)

#monocleで時系列解析
data <- GetAssayData(fib_subset, slot = "counts")
pd <- new("AnnotatedDataFrame", data = data.frame(cell_type = fib_subset@meta.data$original_ann_level_4_refined, row.names = colnames(fib_subset)))
fd <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
fd <- new("AnnotatedDataFrame", data = fd)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "cell_type")

# Pseudotimeプロットの描画（特定の遺伝子発現量を基にしたプロット）
gene_id <- "LRP1"
gene_expression <- as.matrix(cds@assayData$exprs)[gene_id, ]
cds@phenoData@data$GeneExpression <- gene_expression

plot1 <- plot_cell_trajectory(cds, color_by = "GeneExpression", show_branch_points = TRUE, show_tree = TRUE) +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = paste("Expression of", gene_id, "on Pseudotime Trajectory"))

print(plot1)

# 特定の遺伝子発現量を基にしたプロットの描画（ドットサイズを大きくする）
plot2 <- plot_cell_trajectory(cds, color_by = "GeneExpression", show_branch_points = TRUE, show_tree = TRUE, cell_size = 2) +
  scale_color_gradient(low = "lightgrey", high = "darkblue") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(title = paste("Expression of", gene_id, "on Pseudotime Trajectory"))

print(plot2)

#PROGENy解析
DefaultAssay(fib_subset) <- "RNA"  # 正しいアッセイ名に適宜変更
Idents(fib_subset) <- fib_subset@meta.data$reannotated_fibro_type
avg_expr <- AverageExpression(fib_subset, return.seurat = FALSE)$RNA
data("model_full", package = "progeny")
progeny_scores <- progeny(expr = as.matrix(avg_expr), scale = TRUE, organism = "Human", model = model_mouse_full, top = 500)
progeny_scores_t <- t(progeny_scores)
desired_order <- c("Alveolar fibroblasts", "Inflammatory fibroblasts1", "Inflammatory fibroblasts2", "Fibrotic fibroblasts")
my_palette <- colorRampPalette(c("white", "purple"))(100)

# プロット
pheatmap(
  progeny_scores_t[, desired_order],
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = my_palette,
  main = "PROGENy Pathway Activity per Fibroblast Subtype",
  fontsize_row = 10,
  fontsize_col = 10
)


#Connectome解析
# Connectome解析に使用する代表的な細胞タイプを指定
selected_celltypes <- c(
  # マクロファージ
  "Alveolar macrophages",
  
  # 線維芽細胞系（必ず含める）
  "Adventitial fibroblasts",
  "Alveolar fibroblasts",
  "Fibrotic fibroblasts",
  "Inflammatory fibroblasts1",
  "Inflammatory fibroblasts2",
  "Peribronchial fibroblasts",
  
  # 上皮細胞
  "AT2 homeostatic",
  "Transitional AT2",
  "Club",
  "Basal resting",
  "Multiciliated",
  
  # 血管内皮
  "EC aerocyte capillary",
  "EC general capillary",
  
  # T/B細胞など（任意で追加）
  "CD4 T cells",
  "CD8 T cells",
  "B cells"
)

# 上記のcelltypeリストに該当する細胞だけを抽出
data_IP_subset <- subset(data_IP, subset = original_ann_level_4_refined %in% selected_celltypes)
Idents(data_IP_subset) <- "original_ann_level_4_refined"
Idents(data_IP_subset) <- as.factor(data_IP_subset$original_ann_level_4_refined)

# 細胞数が多いクラスターをダウンサンプリング
large_clusters <- c("Alveolar macrophages", "Basal resting", "Multiciliated", "CD4 T cells", "CD8 T cells", "Club", "EC aerocyte capillary", "EC general capillary")
subset_cells <- unlist(lapply(large_clusters, function(celltype) {
  cells <- WhichCells(data_IP_subset, idents = celltype)
  sample(cells, size = floor(length(cells) * 0.5))
}))

# 小さいクラスタはすべて保持
small_cells <- base::setdiff(colnames(data_IP_subset), unlist(lapply(large_clusters, function(ct) WhichCells(data_IP_subset, idents = ct))))

# 統合して subset 作成
reduced_cells <- c(subset_cells, small_cells)
data_IP_reduced_all <- subset(data_IP_subset, cells = reduced_cells)
Idents(data_IP_reduced_all) <- data_IP_reduced_all$original_ann_level_4_refined
data_IP_reduced_all <- NormalizeData(data_IP_reduced_all)
data_IP_reduced_all <- ScaleData(data_IP_reduced_all)
gc()  # メモリ解放
conflicts_prefer(scales::rescale)
conflicts_prefer(clusterProfiler::select)
connectome_subset <- CreateConnectome(
  object = data_IP_reduced_all ,
  species = "human",
  p.values = TRUE,
  min.cells = 10,
  assay = "RNA"
)

connectome_IP_filtered <- FilterConnectome(connectome = connectome_subset, min.pct = 0.8, min.z = 1.0)
CircosPlot(connectome = connectome_IP_filtered, min.pct = 0.2, min.z = 0.05, title = "Connectome from data_IP")

# Lrp1の受容体としての相互作用だけを抽出
conn_Lrp1 <- connectome_subset %>%dplyr::filter(receptor == "LRP1")
CircosPlot(connectome = conn_Lrp1, min.pct = 0.8, min.z = 1.0, title = "Lrp1受容体を介した相互作用")

cells.of.interest <- c(
  "Alveolar fibroblasts", 
  "Inflammatory fibroblasts1", 
  "Inflammatory fibroblasts2", 
  "Fibrotic fibroblasts"
)

CircosPlot(
  connectome = connectome_IP_filtered,
  targets.include = cells.of.interest,
  min.pct = 0.0,
  min.z = 0.0,
  lab.cex = 0.6,
  title = "All interactions toward selected fibroblasts"
)

# フィルタ済connectomeで、各細胞がターゲットになっている件数を確認
table(connectome_IP_filtered$target)
table(data_IP_reduced_all$original_ann_level_4_refined)
conflicted::conflicts_prefer(pheatmap::pheatmap)

CircosPlot(
  connectome = connectome_subset,  # ← FilterConnectomeを使わない元オブジェクト
  targets.include = c(
    "Alveolar fibroblasts", 
    "Inflammatory fibroblasts1", 
    "Inflammatory fibroblasts2", 
    "Fibrotic fibroblasts"
  ),
  min.pct = 0.7,
  min.z = 0.6,
  lab.cex = 0.7,
  title = "All fibroblast types as targets"
)









