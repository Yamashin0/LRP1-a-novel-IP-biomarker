

setwd("C:\\Users\\yamas\\OneDrive\\デスクトップ\\BLM_Seurat")
file_path <- "C:\\Users\\yamas\\OneDrive\\デスクトップ\\BLM_Seurat\\BLM_Seurat.qs"  
data_BLM <- qread(file_path)

DimPlot(data_BLM, group.by = "condition") 
DimPlot(data_BLM, group.by = "orig.ident")  
DimPlot(data_BLM, group.by = "celltype") 
DimPlot(data_BLM, group.by = "celltype_init") 
DimPlot(data_BLM, group.by = "celltype_fine") 

DimPlot(data_BLM, group.by = "celltype", label = TRUE, repel = TRUE) + 
  theme_minimal() + 
  theme(text = element_text(size = 1), 
        legend.position = "none",
        panel.grid.major = element_blank(),  # 大きなグリッドを消す
        panel.grid.minor = element_blank())  # 小さなグリッドを消す

# 特定の遺伝子発現データのプロット
FeaturePlot(data_BLM, features = "Lrp1", min.cutoff = "q10", max.cutoff = "q60", cols = c("lightgrey", "darkred"), pt.size = 1) +
  scale_color_gradientn(colors = c("lightgrey", "#EE3B3B", "#8B2323", "darkred")) +
  theme_minimal() +
  labs(title = "Lrp1 Expression in Fibroblasts",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # 格子を消す
  scale_color_gradientn(colors = c("lightgrey", "#FF4500", "#B22222", "darkred"))  # コントラストを強化

#グレースケール
FeaturePlot(data_BLM, features = "Lrp1", min.cutoff = "q10", max.cutoff = "q60", cols = c("lightgrey", "black"), pt.size = 1) +
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

# 日付をまとめる。
# メタデータを取得
meta_data <- data_BLM@meta.data
meta_data$day_group <- gsub(".*(day[0-9]+).*", "\\1", meta_data$orig.ident)
print(table(meta_data$day_group))
data_BLM <- AddMetaData(object = data_BLM, metadata = meta_data)
head(data_BLM@meta.data)


# Fibroblast細胞のみを取り出してリクラスター。元のメタデータに追加
fibroblasts <- subset(data_BLM, subset = celltype %in% c("Fibroblast", "Fibroblast_Peribronchiolar"))
fibroblasts <- SCTransform(fibroblasts, verbose = FALSE)
fibroblasts <- RunPCA(fibroblasts, verbose = FALSE)
fibroblasts <- RunUMAP(fibroblasts, dims = 1:10, verbose = FALSE)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:10, verbose = FALSE)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.5, verbose = FALSE)
fibroblast_clusters <- Idents(fibroblasts)
data_BLM$fibroblast_clusters <- NA
data_BLM$fibroblast_clusters[Cells(fibroblasts)] <- as.character(fibroblast_clusters)
data_BLM$new_celltype <- ifelse(is.na(data_BLM$fibroblast_clusters), data_BLM$celltype, paste(data_BLM$celltype, data_BLM$fibroblast_clusters, sep = "_"))
data_BLM <- AddMetaData(object = data_BLM, metadata = data_BLM$new_celltype, col.name = "new_celltype")

DimPlot(fibroblasts, pt.size = 1.5) 
DimPlot(fibroblasts, group.by = "condition", pt.size = 1.5) 
DimPlot(fibroblasts, group.by = "day_group", pt.size = 1.5) 

DimPlot(
  fibroblasts, 
  group.by = "fibroblast_clusters", 
  pt.size = 1.5, 
  label = TRUE,       # クラスターラベルを表示
  label.size = 5      # ラベルのサイズを調整
) +
  labs(title = "UMAP with Cluster Labels") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

FeaturePlot(fibroblasts, features = c("Lrp1"), max.cutoff = 3, cols = c("grey", "darkred"), pt.size = 1.5)

FeaturePlot(fibroblasts, features = "Lrp1", min.cutoff = "q10", max.cutoff = "q90", cols = c("lightgrey", "darkred"), pt.size = 1.5) +
  theme_minimal() +
  labs(title = "Lrp1 Expression in Fibroblasts",
       x = "UMAP 1",
       y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # 格子を消す
  scale_color_gradientn(colors = c("lightgrey", "#A8A8A8", "black"), values = c(0, 0.8, 1))  # コントラストを強化


genes_to_plot <- c("Lrp1", "Col1a1", "Col1a2", "Bpifb1", "Sucla2")
dot_plot <- DotPlot(fibroblasts, features = genes_to_plot) +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  theme_minimal() +
  labs(title = "Gene Expression of Lrp1, Col1a1, Col1a2, Bpifb1, Sucla2 by Cluster", 
       x = "Genes", y = "Clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))
print(dot_plot)

# LRP1の発現をfibroblast_clustersごとにバイオリンプロットで表示（SCTアッセイを指定）
VlnPlot(
  fibroblasts, 
  features = "Lrp1",  # LRP1の発現値の遺伝子名
  group.by = "fibroblast_clusters",  # クラスターのカテゴリ
  pt.size = 0,  # ドットのサイズを小さくするか非表示にする
  assay = "SCT"  # SCTアッセイを指定
) +
  scale_fill_manual(values = rep("grey", 14)) +  # すべてのクラスターに同じ色を設定
  labs(title = "LRP1 Expression Across Fibroblast Clusters (SCT Assay)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  # 格子を消す

# Acta2の発現をfibroblast_clustersごとにバイオリンプロットで表示（SCTアッセイを指定）
VlnPlot(
  fibroblasts, 
  features = "Acta2",  # LRP1の発現値の遺伝子名
  group.by = "fibroblast_clusters",  # クラスターのカテゴリ
  pt.size = 0,  # ドットのサイズを小さくするか非表示にする
  assay = "SCT"  # SCTアッセイを指定
) +
  scale_fill_manual(values = rep("grey", 14)) +  # すべてのクラスターに同じ色を設定
  labs(title = "LRP1 Expression Across Fibroblast Clusters (SCT Assay)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  # 格子を消す

# fibroblasts_clustersにおいてLRP1発現データの比較
# LRP1発現データの抽出（filterをdplyrから明示的に呼び出す）
data <- fibroblasts@meta.data %>%
  mutate(Lrp1_expression = fibroblasts@assays$SCT@data["Lrp1", ]) %>%
  dplyr::filter(!is.na(Lrp1_expression))
# Kruskal-Wallis検定の実施
kruskal_test <- kruskal.test(Lrp1_expression ~ fibroblast_clusters, data = data)
# Dunn検定の実施
dunn_test <- dunnTest(Lrp1_expression ~ fibroblast_clusters, data = data, method = "bonferroni")
# Dunn検定の結果をデータフレームに変換
dunn_results <- as.data.frame(dunn_test$res)
# 有意差のあるペアを抽出（selectをdplyrから明示的に呼び出す）
significant_pairs <- dunn_results %>%
  dplyr::filter(P.adj < 0.05) %>%
  dplyr::select(Comparison, P.adj)

# LRP1発現をクラスターごとに可視化（箱ひげ図とドットプロットを重ねて表示）
plot <- ggplot(data, aes(x = fibroblast_clusters, y = Lrp1_expression)) +
  geom_jitter(width = 0.2, alpha = 0.3, color = "grey70") +  # ドットプロットを背景に追加
  geom_boxplot(outlier.size = 0, fill = NA, color = "black") +  # 箱ひげ図を前面に表示
  labs(title = "LRP1 Expression Across Fibroblast Clusters", x = "Clusters", y = "LRP1 Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),  # 主な格子線を消す
    panel.grid.minor = element_blank()   # 小さな格子線を消す
  )

# Dunn検定の結果をチャートに追加
for (i in 1:nrow(significant_pairs)) {
  plot <- plot + geom_signif(
    comparisons = list(strsplit(as.character(significant_pairs$Comparison[i]), " - ")),
    annotations = format(significant_pairs$P.adj[i], digits = 2),
    y_position = max(data$Lrp1_expression) * (1 + 0.05 * i)
  )
}

print(plot)

# Lrp1の発現値をクラスターごとにグループ化し、平均値を計算
cluster_means <- data %>%
  group_by(fibroblast_clusters) %>%
  summarize(mean_Lrp1_expression = mean(Lrp1_expression, na.rm = TRUE))
print(cluster_means)

# クラスター0を基準グループとしてDunnett検定を実行
data <- data %>%
  mutate(fibroblast_clusters = as.factor(fibroblast_clusters))
dunnett_test <- glht(
  aov(Lrp1_expression ~ fibroblast_clusters, data = data),
  linfct = mcp(fibroblast_clusters = "Dunnett")
)
summary(dunnett_test)

# クラスター0と他のすべてのクラスターの比較
markers_cluster_0 <- FindMarkers(
  object = fibroblasts,
  ident.1 = "0",                    # クラスター0を指定
  group.by = "fibroblast_clusters", # クラスタリングの基準
  test.use = "wilcox",              # ウィルコクソン検定を使用
  only.pos = TRUE                   # 発現が上昇している遺伝子のみ
)

# 上昇している遺伝子のトップ300を抽出
top_300_genes <- markers_cluster_0[order(markers_cluster_0$avg_log2FC, decreasing = TRUE), ][1:300, ]
file_path <- "E:/BLM_Seurat/top_300_genes_cluster_0_vs_others.csv"
write.csv(top_300_genes, file = file_path, row.names = TRUE)



# FibroblastとFibroblast_Peribronchiolar細胞のみを取り出しLRP1の発現のフラフを描画
fibroblasts <- subset(data_BLM, subset = celltype %in% c("Fibroblast", "Fibroblast_Peribronchiolar"))
lrp1_data <- FetchData(fibroblasts, vars = c("Lrp1", "day_group", "orig.ident"))
individual_means <- lrp1_data %>%
  group_by(day_group, orig.ident) %>%
  summarise(mean_expression = mean(Lrp1))
summary_stats <- individual_means %>%
  group_by(day_group) %>%
  summarise(mean_expression = mean(mean_expression),
            se_expression = sd(mean_expression) / sqrt(n()))
comparisons <- list(c("day0", "day3"), c("day3", "day7"), c("day7", "day14"), c("day14", "day28"), c("day28", "day42"), c("day42", "day63"))
ggplot(summary_stats, aes(x = day_group, y = mean_expression)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", fill = "transparent") +
  geom_errorbar(aes(ymin = mean_expression - se_expression, ymax = mean_expression + se_expression), 
                width = 0.2, position = position_dodge(0.9)) +
  geom_point(data = individual_means, aes(x = day_group, y = mean_expression), 
             position = position_jitter(width = 0.2), size = 2, color = "black", fill = "white", shape = 21) +
  geom_signif(comparisons = comparisons, map_signif_level = TRUE, test = "t.test") +
  labs(title = "LRP1 Expression by Condition",
       x = "Condition",
       y = "Mean LRP1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))


# FibroblastとFibroblast_Peribronchiolar細胞のみを取り出し各個体のLRP1発現平均値をcsvに保存
fibroblasts <- subset(data_BLM, subset = celltype %in% c("Fibroblast", "Fibroblast_Peribronchiolar"))
lrp1_data <- FetchData(fibroblasts, vars = c("Lrp1", "day_group", "orig.ident"))
individual_means <- lrp1_data %>%
  group_by(day_group, orig.ident) %>%
  summarise(mean_expression = mean(Lrp1))
individual_means_list <- split(individual_means$mean_expression, individual_means$orig.ident)
print(individual_means_list)
write.csv(individual_means_list, "individual_means.csv", row.names = FALSE)

#day28やday42において変動する遺伝子を抽出しcsvファイルに保存する。
Idents(fibroblasts) <- "day_group"
day00_vs_day28 <- FindMarkers(fibroblasts, ident.1 = "day28", ident.2 = "day00")
day00_vs_day42 <- FindMarkers(fibroblasts, ident.1 = "day42", ident.2 = "day00")
top100_genes_day28 <- day00_vs_day28 %>%
  arrange(p_val_adj) %>%
  head(100)
top100_genes_day42 <- day00_vs_day42 %>%
  arrange(p_val_adj) %>%
  head(100)
write.csv(top100_genes_day28, "top100_genes_day28_vs_day00.csv", row.names = TRUE)
write.csv(top100_genes_day42, "top100_genes_day42_vs_day00.csv", row.names = TRUE)

# CytoTRACEの実行
# 発現行列の抽出
expression_matrix <- as.matrix(fibroblasts@assays$RNA@counts)
cytotrace_results <- CytoTRACE(expression_matrix)
plotCytoTRACE(cytotrace_results)
fibroblasts[["CytoTRACE"]] <- cytotrace_results$CytoTRACE
DimPlot(fibroblasts, reduction = "umap", group.by = "CytoTRACE")
FeaturePlot(fibroblasts, features = "CytoTRACE", reduction = "umap", 
            cols = c("lightgrey", "blue"), pt.size = 1.5) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme_minimal() +
  ggtitle("CytoTRACE Score on UMAP")

#cellchat解析
DefaultAssay(data_BLM) <- "RNA"
cellchat <- createCellChat(object = data_BLM, group.by = "combined_labels")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

pathways.show <- c("FGF") # 例としてTGFbパスウェイ
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_heatmap(cellchat)
netVisual_bubble(cellchat, sources.use = 1, targets.use = 2, signaling = pathways.show)
netAnalysis_river(cellchat, signaling = pathways.show)
netAnalysis_contribution(cellchat, signaling = pathways.show)
netAnalysis_pie(cellchat, signaling = pathways.show)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = TRUE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = TRUE, title.name = "Interaction weights/strength")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, slot.name = "netP")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_chord_gene(cellchat, signaling = pathways.show, layout = "chord")

# 使用可能なシグナルパスウェイを取得
pathways.show <- subsetCommunication(cellchat)$pathway_name
print(unique(pathways.show))

netVisual_heatmap(cellchat, signaling = "FGF")
netVisual_heatmap(cellchat, signaling = "PDGF")
netVisual_chord_gene(cellchat, signaling = "FGF")
netVisual_chord_gene(cellchat, signaling = "PDGF")
netVisual_chord_gene(cellchat, signaling = "FGF", title.name = "FGF Signaling Pathway")
netVisual_aggregate(cellchat, signaling = "FGF", layout = "circle")
netVisual_chord_cell(cellchat, signaling = "FGF", title.name = "FGF Signaling Pathway")
netVisual_bubble(cellchat, signaling = "FGF", title.name = "FGF Signaling Pathway")
netVisual_bubble(cellchat, signaling = "PDGF", title.name = "PDGF Signaling Pathway")

# ネットワークのリバウンドプロットを作成
pathways.show <- "COLLAGEN"  
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

pathways.show <- "PDGF"  # 表示するシグナルパスウェイを指定
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#2024年9月11日時系列並べ
str(data_BLM@meta.data)
str(data_BLM@meta.data$day_group)
table(data_BLM@meta.data$celltype)
DimPlot(data_BLM, group.by = "day_group") 
DimPlot(data_BLM, group.by = "celltype_fine") 

# LRP1発現の抽出
# "celltype"が細胞種、"day_group"が時間情報
lrp1_expression <- AverageExpression(data_BLM, features = "Lrp1", group.by = c("celltype", "day_group"))

# データを取り出す
lrp1_data <- as.data.frame(t(lrp1_expression$RNA))

# celltypeとday_groupを '_' で分割
lrp1_data <- as.data.frame(t(lrp1_expression$RNA)) %>%
  rownames_to_column(var = "celltype_day_group") %>%
  separate(celltype_day_group, into = c("celltype", "day_group"), sep = "_", extra = "merge")

# 分離されたデータを確認
head(lrp1_data)

# factorを使用して日付の順番を指定
lrp1_data$day_group <- factor(lrp1_data$day_group, levels = c("day00", "day03", "day07", "day14", "day28", "day42", "day63"))

# dplyrのfilterを使用してNAを削除
lrp1_data <- lrp1_data %>% dplyr::filter(!is.na(day_group))

# ヒートマップの作成（コントラスト強調）
ggplot(lrp1_data, aes(x = day_group, y = celltype, fill = V1)) +
  geom_tile() +
  scale_fill_gradient(low = "lightgrey", high = "darkred", limits = c(0, max(lrp1_data$V1))) +
  labs(title = "Time-series LRP1 Expression in Different Cell Types",
       x = "Time (Days)", 
       y = "Cell Types",
       fill = "LRP1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # X軸ラベルを角度調整

#celltype_fineで区別
# "Fibro_"で始まる細胞種を抽出
# "Fibro_"で始まる細胞種を抽出

# "Fibro_alveolar" を含むすべての細胞種をサブセット化
fibro_alveolar_data <- subset(data_BLM, subset = celltype_fine %in% grep("Fibro_alveolar", data_BLM@meta.data$celltype_fine, value = TRUE))

# "Fibro_alveolar" に統一
fibro_alveolar_data@meta.data$celltype_fine <- "Fibro_alveolar"

# LRP1発現データを抽出
lrp1_fibro_alveolar_expression <- AverageExpression(fibro_alveolar_data, features = "Lrp1", group.by = "day_group")

# データフレームを作成
lrp1_fibro_alveolar_data <- as.data.frame(t(lrp1_fibro_alveolar_expression$RNA))

# データフレームを転置
lrp1_fibro_alveolar_data <- as.data.frame(t(lrp1_fibro_alveolar_expression$RNA))

# 列名にday_groupを適用
colnames(lrp1_fibro_alveolar_data) <- "LRP1_expression"

# day_group列を作成
lrp1_fibro_alveolar_data$day_group <- c("day00", "day03", "day07", "day14", "day28", "day42", "day63")

# celltype_fine列を追加
lrp1_fibro_alveolar_data$celltype_fine <- "Fibro_alveolar"

# ヒートマップの作成
ggplot(lrp1_fibro_alveolar_data, aes(x = day_group, y = celltype_fine, fill = LRP1_expression)) +
  geom_tile() +
  scale_fill_gradient(low = "lightblue", high = "darkred", limits = c(700, 1300), na.value = "lightgrey") +
  labs(title = "Time-series LRP1 Expression in Fibro_alveolar",
       x = "Time (Days)", 
       y = "Fibroblast Subtypes",
       fill = "LRP1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


#2024年9月17日論文用のDEG解析など

# day00とそれ以外の2つのグループに分けて表示.バイオリンプロット
data_BLM@meta.data$day_group_simple <- ifelse(data_BLM@meta.data$day_group == "day00", "day00", "other_days")

# 必要な細胞種のみを選択するリスト
cells_to_keep <- c("Endo_capillary", "Fibroblast", 
                   "Fibroblast_Peribronchiolar", "Epi_AT2", "Epi_AT1", 
                   "Bcell", "MoMac", "Mesothelial", "Epi_ciliated_clara",
                   "SMC", "DC_cDC2",
                   "Pericyte")

# メタデータからcelltype列を使用してフィルタリング
filtered_cells <- data_BLM@meta.data$celltype %in% cells_to_keep
filtered_data_BLM <- subset(data_BLM, cells = which(filtered_cells))

# LRP1発現をcelltype毎に表示、day00とそれ以外の2グループに分けて、格子を消す
VlnPlot(filtered_data_BLM, features = "Lrp1", group.by = "celltype", split.by = "day_group_simple", 
        split.plot = TRUE, pt.size = 0.5) +
  ggtitle("LRP1 Expression Across Cell Types (Day00 vs Other Days)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # 格子を消す
  scale_fill_manual(values = c("day00" = "#FCFCFC", "other_days" = "#828282"))

# LRP1発現をcelltype毎に表示、day00とそれ以外の2グループに分けて、格子を消す
VlnPlot(filtered_data_BLM, features = "Lrp1", group.by = "celltype", split.by = "day_group_simple", 
        split.plot = TRUE, pt.size = 0) +  # pt.size = 0 でドットを非表示
  ggtitle("LRP1 Expression Across Cell Types (Day00 vs Other Days)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # 格子を消す
  scale_fill_manual(values = c("day00" = "#FCFCFC", "other_days" = "#828282"))

#細胞種毎にLRP1の発現をdot plotで示す
# 必要な細胞種のみを選択するリスト
# 細胞の順番を指定する
cell_order <- c("Bcell", "DC_cDC2", "Endo_capillary", "Epi_AT1", "Epi_AT2", 
                "Epi_ciliated_clara", "Fibroblast", "Fibroblast_Peribronchiolar", 
                "Mesothelial", "MoMac", "Pericyte", "SMC")

# メタデータからcelltype列を使用してフィルタリング
filtered_cells <- data_BLM@meta.data$celltype %in% cell_order
filtered_data_BLM <- subset(data_BLM, cells = which(filtered_cells))

# LRP1発現が0以上の細胞のみを表示
# Dot PlotでLRP1の発現を細胞種毎に表示、色を指定し、格子を消してドットを大きく
DotPlot(filtered_data_BLM, features = "Lrp1", group.by = "celltype", cols = c("#FCFCFC", "black"), dot.scale = 10) +
  ggtitle("LRP1 Expression Across Selected Cell Types") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # X軸のラベルを45度回転
        panel.grid.major = element_blank(),  # 格子（大きなグリッド）を消す
        panel.grid.minor = element_blank())  # 格子（小さなグリッド）を消す


#細胞種毎にHCvsBLMでLRP1発現の有意差検定
# 正しいLrp1の行名を指定
lrp1_name <- "Lrp1"  # "Lrp1"が正しい遺伝子名

# 正しいLrp1の行名を指定
lrp1_name <- "Lrp1"  # "Lrp1"が正しい遺伝子名

# 細胞種ごとにday00とother_daysの発現データを抽出し、Wilcoxon rank-sum testを実行
celltypes <- unique(filtered_data_BLM@meta.data$celltype)
p_values <- data.frame(celltype = character(), p_value = numeric(), stringsAsFactors = FALSE)

for (celltype in celltypes) {
  # day00とother_daysのセルをmeta.dataから取得
  day00_cells <- rownames(filtered_data_BLM@meta.data[filtered_data_BLM@meta.data$day_group_simple == "day00" & 
                                                        filtered_data_BLM@meta.data$celltype == celltype, ])
  
  other_days_cells <- rownames(filtered_data_BLM@meta.data[filtered_data_BLM@meta.data$day_group_simple == "other_days" & 
                                                             filtered_data_BLM@meta.data$celltype == celltype, ])
  
  # データが存在するか確認して統計検定を実行
  if (length(day00_cells) > 0 & length(other_days_cells) > 0) {
    day00_data <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = day00_cells)
    other_days_data <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = other_days_cells)
    
    # Wilcoxon rank-sum testを実行
    p_value <- wilcox.test(day00_data[[lrp1_name]], other_days_data[[lrp1_name]])$p.value
    
    # p-valueをデータフレームに保存
    p_values <- rbind(p_values, data.frame(celltype = celltype, p_value = p_value))
  }
}

# p-valueの結果を表示
p_values

# 各細胞種ごとにday00とother_daysの平均値を計算
mean_values <- data.frame(celltype = character(), mean_day00 = numeric(), mean_other_days = numeric(), stringsAsFactors = FALSE)

# 有意差があった細胞種
significant_celltypes <- c("Bcell", "Fibroblast", "SMC", "Neutrophil", "Epi_ciliated_clara", "Pericyte")

for (celltype in significant_celltypes) {
  # day00とother_daysのセルをmeta.dataから取得
  day00_cells <- rownames(filtered_data_BLM@meta.data[filtered_data_BLM@meta.data$day_group_simple == "day00" & 
                                                        filtered_data_BLM@meta.data$celltype == celltype, ])
  
  other_days_cells <- rownames(filtered_data_BLM@meta.data[filtered_data_BLM@meta.data$day_group_simple == "other_days" & 
                                                             filtered_data_BLM@meta.data$celltype == celltype, ])
  
  # 平均値を計算
  if (length(day00_cells) > 0 & length(other_days_cells) > 0) {
    day00_data <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = day00_cells)
    other_days_data <- FetchData(filtered_data_BLM, vars = lrp1_name, cells = other_days_cells)
    
    # 平均値を計算
    mean_day00 <- mean(day00_data[[lrp1_name]])
    mean_other_days <- mean(other_days_data[[lrp1_name]])
    
    # 平均値をデータフレームに保存
    mean_values <- rbind(mean_values, data.frame(celltype = celltype, mean_day00 = mean_day00, mean_other_days = mean_other_days))
  }
}

# 結果を表示
mean_values


#2024年9月19日
# PROGENyを使用してscRNA-seqデータを解析
# PROGENy解析の実行・SCT正規化されたデータを使用して解析を実行・SCTデータを疎行列から標準的な行列形式に変換
dense_matrix <- as.matrix(fibroblasts@assays$SCT@data)

progeny_scores <- progeny(dense_matrix, 
                          scale = TRUE,  # データのスケーリング
                          organism = "Mouse", 
                          top = 100,     # 使用する上位ターゲット遺伝子数
                          perm = 1000)   # パーミュテーションの回数

Heatmap(progeny_scores)

# progeny_scoresの列名に対応する細胞タイプを取得
celltype_annotation <- fibroblasts@meta.data[colnames(progeny_scores), "celltype"]

# 注釈付きヒートマップの作成# クラスタリング
Heatmap(progeny_scores, 
        show_row_names = FALSE, 
        show_column_names = TRUE,
        row_title = "Pathways",
        column_title = "Cells",
        top_annotation = HeatmapAnnotation(
          celltype = celltype_annotation,  # 修正された注釈データ
          col = list(celltype = c("Fibroblast" = "blue", "Myofibroblast" = "red"))
        ))

progeny_clusters <- hclust(dist(t(progeny_scores)))
plot(progeny_clusters)

# パスウェイのスコアを抽出・UMAPでの描出・転置は必要なければ行わない。
#転置
progeny_scores <- t(progeny_scores)
rownames(progeny_scores)

#各スコアをfibroblastsのメタデータに記入
mapk_scores <- progeny_scores["MAPK", ]
fibroblasts$MAPK_score <- mapk_scores

WNT_scores <- progeny_scores["WNT", ]
fibroblasts$WNT_score <- WNT_scores

Androgen_scores <- progeny_scores["Androgen", ]
fibroblasts$ANdrogen_score <- Androgen_scores

EGFR_scores <- progeny_scores["EGFR", ]
fibroblasts$EGFR_score <- EGFR_scores

Estrogen_scores <- progeny_scores["Estrogen", ]
fibroblasts$Estrogen_score <- Estrogen_scores

Hypoxia_scores <- progeny_scores["Hypoxia", ]
fibroblasts$Hypoxia_score <- Hypoxia_scores

JAK.STAT_scores <- progeny_scores["JAK.STAT", ]
fibroblasts$JAK.STAT_score <- JAK.STAT_scores

NFkB_scores <- progeny_scores["NFkB", ]
fibroblasts$NFkB_score <- NFkB_scores

p53_scores <- progeny_scores["p53", ]
fibroblasts$p53_score <- p53_scores

PI3K_scores <- progeny_scores["PI3K", ]
fibroblasts$PI3K_score <- PI3K_scores

TGFb_scores <- progeny_scores["TGFb", ]
fibroblasts$TGFb_score <- TGFb_scores

TNFa_scores <- progeny_scores["TNFa", ]
fibroblasts$TNFa_score <- TNFa_scores

Trail_scores <- progeny_scores["Trail", ]
fibroblasts$Trail_score <- Trail_scores

VEGF_scores <- progeny_scores["VEGF", ]
fibroblasts$VEGF_score <- VEGF_scores


# MAPKパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
umap_results <- umap::umap(t(progeny_scores))
plot(umap_results$layout, col = rainbow(length(unique(fibroblasts@meta.data$celltype)))[factor(fibroblasts@meta.data$celltype)])

high_mapk <- colnames(progeny_scores)[mapk_scores > quantile(mapk_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "MAPK_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of MAPK Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# WNTパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
high_WNT <- colnames(progeny_scores)[WNT_scores > quantile(WNT_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "WNT_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of WNT Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# PI3Kパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
high_PI3K <- colnames(progeny_scores)[PI3K_scores > quantile(PI3K_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "PI3K_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of PI3K Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# EGFRパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
high_EGFR <- colnames(progeny_scores)[EGFR_scores > quantile(EGFR_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "EGFR_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of EGFR Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# JAK.STATパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
high_JAK.STAT <- colnames(progeny_scores)[JAK.STAT_scores > quantile(JAK.STAT_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "JAK.STAT_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of JAK.STAT Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Estrogenパスウェイのスコアを抽出高スコアの細胞を抽出（例: 上位25%の細胞）高スコア細胞の可視化
high_Estrogen <- colnames(progeny_scores)[Estrogen_scores > quantile(Estrogen_scores, 0.9)]

FeaturePlot(
  fibroblasts,
  features = "Estrogen_score",
  reduction = "umap",
  cols = c("lightgrey", "darkred"),  # カラースケールを調整
  min.cutoff = "q10",  # 下限のカットオフを設定
  max.cutoff = "q90",  # 上限のカットオフを設定
  pt.size = 1  # プロットサイズ
) +
  ggtitle("UMAP of Estrogen Pathway Scores") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

unique(fibroblasts@meta.data$day_group)
unique(fibroblasts@meta.data$fibroblast_clusters)

# day00とday28のデータを抽出
day00_data <- subset(fibroblasts, subset = day_group == "day00")
day28_data <- subset(fibroblasts, subset = day_group == "day28")

# WNTスコアの比較
wnt_scores_day00 <- day00_data@meta.data$WNT_score
wnt_scores_day28 <- day28_data@meta.data$WNT_score

# MAPKスコアの比較
mapk_scores_day00 <- day00_data@meta.data$MAPK_score
mapk_scores_day28 <- day28_data@meta.data$MAPK_score

# PI3Kスコアの比較
pi3k_scores_day00 <- day00_data@meta.data$PI3K_score
pi3k_scores_day28 <- day28_data@meta.data$PI3K_score

# WNTスコアのt検定
wnt_test <- t.test(wnt_scores_day00, wnt_scores_day28)
print(wnt_test)

# MAPKスコアのt検定
mapk_test <- t.test(mapk_scores_day00, mapk_scores_day28)
print(mapk_test)

# PI3Kスコアのt検定
pi3k_test <- t.test(pi3k_scores_day00, pi3k_scores_day28)
print(pi3k_test)

# WNTスコアのボックスプロット
ggplot(fibroblasts@meta.data, aes(x = day_group, y = WNT_score)) +
  geom_boxplot(aes(fill = day_group)) +
  labs(title = "WNT Score Comparison", x = "Day Group", y = "WNT Score") +
  theme_minimal()

# MAPKスコアのボックスプロット
ggplot(fibroblasts@meta.data, aes(x = day_group, y = MAPK_score)) +
  geom_boxplot(aes(fill = day_group)) +
  labs(title = "MAPK Score Comparison", x = "Day Group", y = "MAPK Score") +
  theme_minimal()

# PI3Kスコアのボックスプロット
ggplot(fibroblasts@meta.data, aes(x = day_group, y = PI3K_score)) +
  geom_boxplot(aes(fill = day_group)) +
  labs(title = "PI3K Score Comparison", x = "Day Group", y = "PI3K Score") +
  theme_minimal()

# day_groupがday00とday28のみを抽出
filtered_data <- fibroblasts@meta.data %>%
  dplyr::filter(day_group %in% c("day00", "day28"))

# PI3Kスコアのボックスプロットとt検定結果の表示
ggplot(filtered_data, aes(x = day_group, y = PI3K_score)) +
  geom_boxplot(aes(fill = day_group)) +
  labs(title = "PI3K Score Comparison", x = "Day Group", y = "PI3K Score") +
  theme_minimal() +
  stat_compare_means(method = "t.test", label = "p.signif")  # 有意差検定の結果を追加

# LRP1発現値をmeta.dataに追加
fibroblasts@meta.data$LRP1_expression <- fibroblasts@assays$SCT@data["Lrp1", ]

# LRP1の発現値に基づいて細胞を「高い」「低い」に分ける（中央値で分ける例）
fibroblasts@meta.data$LRP1_group <- ifelse(fibroblasts@meta.data$LRP1_expression < median(fibroblasts@meta.data$LRP1_expression), "Low", "High")

# 正しいパスウェイスコアの名前をリスト化
pathways <- c("ANdrogen_score", "EGFR_score", "Estrogen_score", "Hypoxia_score", "JAK.STAT_score", 
              "MAPK_score", "NFkB_score", "p53_score", "PI3K_score", "TGFb_score", 
              "TNFa_score", "Trail_score", "VEGF_score", "WNT_score") 

# DotPlotの作成
DotPlot(
  fibroblasts,
  features = pathways,
  group.by = "LRP1_group"
) + 
  theme_minimal() + 
  labs(title = "Pathway Scores by LRP1 Expression Group") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # X軸のラベルを45度回転
    panel.grid.major = element_blank(),  # 格子を消す
    panel.grid.minor = element_blank()
  ) +
  scale_color_gradient(low = "grey", high = "darkred") +  # 色のカラースケールを変更
  scale_size(range = c(3, 8))  # ドットのサイズ範囲を変更してコントラストをはっきり

#箱ひげ図
plots <- lapply(pathways, function(pathway) {
  ggplot(fibroblasts@meta.data, aes_string(x = "LRP1_group", y = pathway)) +
    geom_boxplot(aes(fill = LRP1_group)) +
    labs(title = paste(pathway, "Score Comparison by LRP1 Expression"), x = "LRP1 Group", y = paste(pathway, "Score")) +
    theme_minimal() +
    stat_compare_means(method = "t.test", label = "p.signif")  # 有意差検定の結果を表示
})

# プロットの表示
print(plots)

# Wilcoxon検定を例に、各パスウェイについてp値を計算
# p値を格納するデータフレームを作成
p_values <- data.frame(Pathway = pathways, p_value = NA)

# データの抽出確認
# LRP1_groupのラベルを統一
fibroblasts@meta.data$LRP1_group <- tolower(fibroblasts@meta.data$LRP1_group)

# p値を格納するデータフレームを作成
p_values <- data.frame(Pathway = pathways, p_value = NA)

# LRP1のグループに基づいて各パスウェイのスコアを比較
for (pathway in pathways) {
  high_group <- fibroblasts@meta.data[fibroblasts@meta.data$LRP1_group == "high", pathway, drop = FALSE]
  low_group <- fibroblasts@meta.data[fibroblasts@meta.data$LRP1_group == "low", pathway, drop = FALSE]
  
  # NAを除外
  high_group <- high_group[!is.na(high_group)]
  low_group <- low_group[!is.na(low_group)]
  
  # 十分なデータがあるかをチェックしてからテストを実行
  if (length(high_group) > 1 & length(low_group) > 1) {
    p_values[p_values$Pathway == pathway, "p_value"] <- wilcox.test(high_group, low_group)$p.value
  } else {
    p_values[p_values$Pathway == pathway, "p_value"] <- NA
  }
}

# p値の確認
print(p_values)

# DotPlotの作成
ggplot(p_values, aes(y = reorder(Pathway, -p_adjust), x = -log10(p_adjust), size = -log10(p_value), color = p_adjust)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +  # p値に基づく色分け
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # 格子を消す
  labs(
    title = "Progeny with Multiple Signaling Pathways",
    y = "Pathway",
    x = "-log10(Adjusted P-value)"
  ) +
  scale_size(range = c(2, 10))  # Dotのサイズ範囲を調整

#fibroblastsの中でLrp1が最も高いクラスター0と最も低いクラスター11を比較する
# クラスター0とクラスター11のデータを抽出
cluster_0_data <- fibroblasts@meta.data %>%
  dplyr::filter(fibroblast_clusters == "0")
cluster_11_data <- fibroblasts@meta.data %>%
  dplyr::filter(fibroblast_clusters == "11")

# 正しいパスウェイスコアの名前をリスト化
pathways <- c("ANdrogen_score", "EGFR_score", "Estrogen_score", "Hypoxia_score", "JAK.STAT_score", 
              "MAPK_score", "NFkB_score", "p53_score", "PI3K_score", "TGFb_score", 
              "TNFa_score", "Trail_score", "VEGF_score", "WNT_score") 

# パスウェイごとに比較し、p値を取得
comparison_results <- data.frame(Pathway = pathways, avg_logFC = NA, p_value = NA)

for (i in seq_along(pathways)) {
  pathway <- pathways[i]
  
# クラスター間でのWilcoxon検定
test_result <- wilcox.test(cluster_0_data[[pathway]], cluster_11_data[[pathway]])
  
# 平均のフォールド変化とp値の記録
  avg_logFC <- mean(cluster_0_data[[pathway]]) - mean(cluster_11_data[[pathway]])
  comparison_results$avg_logFC[i] <- avg_logFC
  comparison_results$p_value[i] <- test_result$p.value
}

# p値の調整 (FDR)
comparison_results$p_adjust <- p.adjust(comparison_results$p_value, method = "fdr")

# 有意な結果をフィルタリング（例：p値 < 0.05）
significant_results <- comparison_results %>% dplyr::filter(p_adjust < 0.05)

# p値で並べ替え
significant_results <- significant_results %>% arrange(p_adjust)

# データフレームを転置
transposed_results <- t(significant_results)

# データフレームとして再変換し、行名を列に変換
transposed_results <- as.data.frame(transposed_results)
colnames(transposed_results) <- transposed_results[1, ]  # 1行目を列名に
transposed_results <- transposed_results[-1, ]  # 1行目を削除

# 行名を列に追加して確認
transposed_results$Pathway <- rownames(transposed_results)
rownames(transposed_results) <- NULL

# Dot Plotの作成（修正版）
ggplot(significant_results, aes(x = reorder(Pathway, p_adjust), y = -log10(p_adjust))) +
  geom_point(aes(size = abs(avg_logFC), color = -log10(p_adjust))) +  # avg_logFCを使用
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Significant Pathways in Cluster 0 vs Cluster 11",
       x = "Pathway",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_size_continuous(range = c(2, 8))  # 点の大きさを調整


# scMetabolismを用いた代謝経路の活性スコアの計算
# メモリ整頓
rm(filtered_data_BLM, individual_means)
gc()
options(future.globals.maxSize = 20 * 1024^3)  # 20GBに設定


# igraphパッケージのdegree関数を優先する
# Ensemblからマウスとヒトの遺伝子の変換を行う
options(timeout = 300)
ensembl <- useEnsembl(biomart = "ensembl")

ensembl <- useEnsembl(biomart = "ensembl", mirror = "www") 
mouse <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# マウス遺伝子をヒト遺伝子に変換
gene_conversion <- getLDS(attributes = c("mgi_symbol"), 
                          filters = "mgi_symbol", 
                          values = mouse_genes, 
                          mart = mouse, 
                          attributesL = c("hgnc_symbol"), 
                          martL = human, 
                          uniqueRows = TRUE)





# ヒトの遺伝子名に変換
rownames(filtered_fibroblasts@assays$RNA@data) <- gene_conversion$HGNC.symbol[match(rownames(filtered_fibroblasts@assays$RNA@data), gene_conversion$MGI.symbol)]

conflicts_prefer(igraph::degree)
small_fibroblasts <- subset(fibroblasts, cells = sample(colnames(fibroblasts), 500))  # 500細胞をランダムに抽出
filtered_fibroblasts <- subset(small_fibroblasts, features = rownames(small_fibroblasts)[1:1000])
countexp.Seurat<-sc.metabolism.Seurat(obj = filtered_fibroblasts, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

# フィルタリングして遺伝子数を減らす
# sc.metabolism.Seuratの実行

small_fibroblasts <- subset(fibroblasts, cells = sample(colnames(fibroblasts), 500))  # 500細胞をランダムに抽出
# 事前に遺伝子数をフィルタリング
filtered_fibroblasts <- subset(small_fibroblasts, features = rownames(small_fibroblasts)[1:1000])

# sc.metabolism.Seuratの実行
fibroblasts_metabolism <- sc.metabolism.Seurat(filtered_fibroblasts, method = "VISION")

fibroblasts_metabolism <- sc.metabolism.Seurat(small_fibroblasts, method = "VISION")

filtered_fibroblasts <- subset(filtered_fibroblasts, features = rownames(filtered_fibroblasts)[1:1000])
fibroblasts_metabolism <- sc.metabolism.Seurat(filtered_fibroblasts, method = "VISION")


# igraphパッケージのdegree関数を優先して使用（どちらか）
conflicts_prefer(igraph::degree)
fibroblasts_metabolism<-sc.metabolism.Seurat(obj = fibroblasts, method = "AUCell", imputation = F, ncores = 2, metabolism.type = "KEGG")

# Glycolysis / Gluconeogenesis のUMAPプロット
p1 <- DimPlot.metabolism(obj = fibroblasts_metabolism, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1.5) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
print(p1)

# Citrate cycle (TCA cycle) のUMAPプロット
p2 <- DimPlot.metabolism(obj = fibroblasts_metabolism, pathway = "Citrate cycle (TCA cycle)", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1.5) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
print(p2)

# Oxidative phosphorylation のUMAPプロット
p3 <- DimPlot.metabolism(obj = fibroblasts_metabolism, pathway = "Oxidative phosphorylation", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1.5) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
print(p3)

# Oxidative phosphorylation のUMAPプロット
p4 <- DimPlot.metabolism(obj = fibroblasts_metabolism, pathway = "Nitrogen metabolism", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1.5) +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
print(p4)

input.pathway<-c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(obj = fibroblasts_metabolism, pathway = input.pathway, phenotype = "seurat_clusters", norm = "y")
BoxPlot.metabolism(obj = fibroblasts_metabolism, pathway = input.pathway, phenotype = "original_ann_level_4", ncol = 1)

# パスウェイリストの拡張
input.pathway <- c(
  "Glycolysis / Gluconeogenesis", 
  "Oxidative phosphorylation", 
  "Citrate cycle (TCA cycle)",
  "Pentose phosphate pathway",           # 追加のパスウェイ例
  "Fatty acid metabolism",               # 追加のパスウェイ例
  "Purine metabolism",                   # 追加のパスウェイ例
  "Pyrimidine metabolism",               # 追加のパスウェイ例
  "Steroid biosynthesis"                 # 追加のパスウェイ例
)

# DotPlotの作成
DotPlot.metabolism(
  obj = fibroblasts_metabolism, 
  pathway = input.pathway, 
  phenotype = "seurat_clusters", 
  norm = "y"
)


# 複数の遺伝子発現をDotPlotで可視化
genes <- c("LRP1", "COL15A1", "CXCL14", "DIO2", "MFAP5", "SCARA5", "PI16", "SPINT2", 
           "LIMCH1", "FGFR4", "MARCKSL1", "APOC1", "MMP23B", "LMOD1", "ATP1B1", "TYRP1", "ACTA2")

# 表示したい遺伝子リストを作成
genes <- c("LRP1", "COL1A2", "DCN", "MFAP4", "LUM", "COL6A3", "CFD", "LMOD1", "ATP1B1", "TYRP1", "ACTA2", "TAGLN", "MYH11", "PDGFRB", "IL1B", "TGFB1", "FN1", "MMP2", "ELN", "SERPINE1")

# HLCA線維芽細胞（Stroma:COL1A2 DCN MFAP4 Fibroblast lineage:LUM COL6A3 CFD Peribronchial fibroblasts:COL15A1 CXCL14 DIO2 
#Adventitial fibroblasts:MFAP5 SCARA5 PI16 Alveolar fibroblasts:SPINT2 LIMCH1 FGFR4  subpleural fibroblasts:MARCKSL1 APOC1 MMP23B
#Myofibroblasts:LMOD1 ATP1B1 TYRP1 
genes <- c("LRP1", "COL1A2", "DCN", "MFAP4", "LUM", "COL6A3", "CFD", "COL15A1", "CXCL14", "DIO2", "MFAP5", "SCARA5", "PI16", "SPINT2", 
           "LIMCH1", "FGFR4", "MARCKSL1", "APOC1", "MMP23B", "LMOD1", "ATP1B1", "TYRP1", "ACTA2")

# 複数の色を使ったカラースケールをカスタマイズしたDotPlot
DotPlot(fibroblasts_data, features = genes, group.by = "seurat_clusters") +
  scale_color_gradientn(colors = c("#CDBA96", "#FFE1FF", "#8B2252")) +  # 青から緑、黄色、赤へのグラデーション
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#2024年9月9日hdWGCNA解析
fibroblasts <- SetupForWGCNA(
  seurat_obj = fibroblasts,
  gene_select = "fraction", # 遺伝子選択方法
  fraction = 0.05, # 発現する必要のある細胞の割合
  wgcna_name = "IP_analysis" # hdWGCNA実験の名前
)

fibroblasts <- MetacellsByGroups(
  seurat_obj = fibroblasts,
  group.by = c("fibroblast_clusters"), # グループ化の基準
  reduction = 'pca', # 利用可能な次元削減
  k = 25, # k近傍パラメータ
  max_shared = 10, # メタセル間で共有できる最大セル数
  ident.group = 'fibroblast_clusters' # メタセルSeuratオブジェクトのIdents設定
)

fibroblasts <- NormalizeMetacells(fibroblasts)

fibroblasts <- SetDatExpr(
  seurat_obj = fibroblasts,
  group_name = "0", # 関心のあるグループ
  group.by = 'fibroblast_clusters', # メタデータ列
  assay = 'SCT', # 使用するアッセイ
  slot = 'data' # 正規化されたデータ
)

fibroblasts <- TestSoftPowers(
  seurat_obj = fibroblasts,
  networkType = 'signed' # ネットワークタイプ
)

# 結果をプロット
plot_list <- PlotSoftPowers(fibroblasts)
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(fibroblasts)
head(power_table)

# construct co-expression network:
fibroblasts <- ConstructNetwork(
  fibroblasts,
  tom_name = 'Fibroblasts', # トポロジカルオーバーラップマトリックスの名前
  overwrite_tom = TRUE # TOMを上書きする
)

PlotDendrogram(fibroblasts, main='Fibroblasts hdWGCNA Dendrogram')
TOM <- GetTOM(fibroblasts)

# Seuratオブジェクトから正規化されたデータを抽出
expr_data <- GetAssayData(fibroblasts, slot = "data")

# 行が遺伝子、列がサンプルになるように転置
datExpr <- t(as.matrix(expr_data))

# ソフトパワーをテスト
powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 結果をプロット
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     type="n", xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit", 
     main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, col="red")

# 共発現ネットワークの構築とモジュール検出
# cor関数の競合を解消する
conflicts_prefer(WGCNA::cor)

# 共発現ネットワークの構築とモジュール検出
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, saveTOMFileBase = "TOM",
                       verbose = 3)

# モジュールカラーの取得
module_colors <- labels2colors(net$colors)

# LRP1が属するモジュールを確認
lrp1_module <- module_colors["LRP1"]
lrp1_module


# LRP1が属するモジュール内の遺伝子を取得
genes_in_lrp1_module <- names(module_colors[module_colors == lrp1_module])
head(genes_in_lrp1_module)

# LRP1と共発現する遺伝子のヒートマップを作成
DoHeatmap(fibroblasts, features = genes_in_lrp1_module[1:10], group.by = "fibroblast_clusters")


#2024年9月26日コネクトーム解析
# Connectomeオブジェクトの作成（種を指定）
species <- "mouse"  # 使用する種を指定
conflicts_prefer(scales::rescale)
connectome_obj <- CreateConnectome(data_BLM, clusters = yamamoto_20240926, species = species)

# Connectomeオブジェクトが正しく作成されたか確認
str(connectome_obj)

# Change density plot line colors by groups
p1 <- ggplot(connectome_obj, aes(x=ligand.scale)) + geom_density() + ggtitle('Ligand.scale')
p2 <- ggplot(connectome_obj, aes(x=recept.scale)) + geom_density() + ggtitle('Recept.scale')
p3 <- ggplot(connectome_obj, aes(x=percent.target)) + geom_density() + ggtitle('Percent.target')
p4 <- ggplot(connectome_obj, aes(x=percent.source)) + geom_density() + ggtitle('Percent.source')
plot_grid(p1,p2,p3,p4)

p1 <- NetworkPlot(connectome_obj,features = 'VEGFA',min.pct = 0.0001,weight.attribute = 'weight_sc',include.all.nodes = T)
p2 <- NetworkPlot(connectome_obj,features = 'VEGFA',min.pct = 0.0001,weight.attribute = 'weight_sc',include.all.nodes = T)

plot_grid(p1,p2,nrow=1)


conflicts_prefer(cowplot::get_legend)
Centrality(connectome_obj,
           modes.include = NULL,
           min.z = 0,  # min.z = 0 を推奨されているため追加
           weight.attribute = 'weight_sc',
           group.by = 'mode')


Centrality(connectome_obj,
           modes.include = c('VEGF'),
           weight.attribute = 'weight_sc',
           min.z = 0,
           group.by = 'mechanism')

p1 <- CellCellScatter(connectome_obj,sources.include = 'Fibroblast',targets.include = 'MoMac',
                      label.threshold = 3,
                      weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 2)
p1 <- p1 + xlim(0,NA) + ylim(0,NA)
p1

p2 <- SignalScatter(connectome_obj, features = c('JAG1','JAG2','LRP1'),label.threshold = 1,weight.attribute = 'weight_sc',min.z = 2)
p2 <- p2 + xlim(2,NA) + ylim(2,NA)
p2

# 名前の衝突を解決
library(conflicted)
conflict_prefer("select", "dplyr")
conflicts_prefer(ComplexHeatmap::draw)

test <- connectome_obj
test <- data.frame(test %>% group_by(vector) %>% top_n(5,weight_sc))


# Connectome解析の実行

#fibroblast,MoMac,Endo_capillaryに注目
cells.of.interest <- c('Fibroblast','MoMac','Endo_capillary')
circos.par(gap.degree = 2)

# Using edgeweights from normalized slot:
# min.pct の設定を緩めて、フィルタリング後にエッジが残るようにする
CircosPlot(
  test,
  weight.attribute = 'weight_norm',
  sources.include = cells.of.interest,
  targets.include = cells.of.interest,
  lab.cex = 1.4,
  min.pct = 0.1,  # 例: フィルタリング基準を緩和
  title = 'Edgeweights from normalized slot'
)


#fibroblasts,MoMac,Epi_ciliated_clara,Epi_AT2に注目
cells.of.interest <- c('Fibroblast','MoMac','Epi_ciliated_clara','Epi_AT2')
circos.par(gap.degree = 2)

# Using edgeweights from normalized slot:
# min.pct の設定を緩めて、フィルタリング後にエッジが残るようにする
CircosPlot(
  test,
  weight.attribute = 'weight_norm',
  sources.include = cells.of.interest,
  targets.include = cells.of.interest,
  lab.cex = 1.2,
  min.pct = 0.1,  # 例: フィルタリング基準を緩和
  title = 'Edgeweights from normalized slot'
)

CircosPlot(test,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from normalized slot')

# Using edgeweights from scaled slot:
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from scaled slot')
# Using separate ligand and receptor expression (from normalized slot)
CircosPlot(test,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')

# Using separate ligand and receptor expression (from scaled slot)
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')

CircosPlot(test,targets.include = 'Fibroblast',lab.cex = 0.6)

CircosPlot(test,sources.include = 'Fibroblast',lab.cex = 0.6)

cells.of.interest <- c('Fibroblast','Epi_AT2')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 1.4,title = 'Edgeweights from scaled slot')

cells.of.interest <- c('Fibroblast','MoMac')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 1.4,title = 'Edgeweights from scaled slot')

cells.of.interest <- c('Fibroblast','Endo_capillary')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 1.2,title = 'Edgeweights from scaled slot')

cells.of.interest <- c('Fibroblast','Epi_AT1')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 1.0,title = 'Edgeweights from scaled slot')

cells.of.interest <- c('Fibroblast','Epi_ciliated_clara')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 1.0,title = 'Edgeweights from scaled slot')

cells.of.interest <- c('Fibroblasts','Submucosal Secretory')
CircosPlot(test,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from scaled slot')










#GSEA解析スクリプト

# Fibroblasts（線維芽細胞）のサブセット化
fibroblast_data <- subset(data_IP_sampled, subset = singlr_labels == "Fibroblasts")

# normallung（健常肺）データのサブセット化
normallung_data <- subset(fibroblast_data, subset = sample_type == "normallung")

# IPFlung（IPF肺）データのサブセット化
IPFlung_data <- subset(fibroblast_data, subset = sample_type == "IPFlung")

# 線維芽細胞データをnormal lungとIPF lungに分ける
normallung_data <- subset(fibroblast_data, subset = sample_type == "normallung")
IPFlung_data <- subset(fibroblast_data, subset = sample_type == "IPFlung")

# LRP1の発現に基づいてグループ分け
normallung_data$LRP1_group <- ifelse(normallung_data@assays$SCT@data["LRP1", ] > median(normallung_data@assays$SCT@data["LRP1", ]), "LRP1_high", "LRP1_low")

# IPF lungのデータにも同様にグループ分け
IPFlung_data$LRP1_group <- ifelse(IPFlung_data@assays$SCT@data["LRP1", ] > median(IPFlung_data@assays$SCT@data["LRP1", ]), "LRP1_high", "LRP1_low")

# Normal lungのグループをIdentsに設定
Idents(normallung_data) <- "LRP1_group"

# IPF lungのグループもIdentsに設定
Idents(IPFlung_data) <- "LRP1_group"

# IPF_LRP1_highとIPF_LRP1_lowに基づくIDの設定
Idents(fibroblast_data) <- "sample_type"

# まずは、IPF_LRP1_highとしてIPFlung内のLRP1_highをサブセット化し、それを一つのグループとします。
fibroblast_data$IPF_LRP1_high <- ifelse(fibroblast_data@assays$SCT@data["LRP1", ] > median(IPFlung_data@assays$SCT@data["LRP1", ]) & fibroblast_data$sample_type == "IPFlung", "IPF_LRP1_high", "other")

# 健常肺を全体として設定
fibroblast_data$IPF_LRP1_high <- ifelse(fibroblast_data$sample_type == "normallung", "normallung", fibroblast_data$IPF_LRP1_high)

# Identsに新しいグループを設定
Idents(fibroblast_data) <- "IPF_LRP1_high"

# IPF_LRP1_highと健常肺全体を比較
DEGs_comparison <- FindMarkers(fibroblast_data, ident.1 = "IPF_LRP1_high", ident.2 = "normallung", logfc.threshold = 0.25)

# GSEA用のランキングリストの作成
ranked_genes_comparison <- DEGs_comparison$avg_log2FC
names(ranked_genes_comparison) <- rownames(DEGs_comparison)

# 降順にソート
ranked_genes_comparison <- sort(ranked_genes_comparison, decreasing = TRUE)

# GSEAの実行
gsea_result_comparison <- GSEA(ranked_genes_comparison, TERM2GENE = pathway_db, pvalueCutoff = 0.05)

# GSEA結果の確認
head(gsea_result_comparison@result)

# dotplotでGSEA結果を可視化
dotplot(gsea_result_comparison) + ggtitle("GSEA of IPF LRP1_high vs Normal Lung")

# enrichplotでGSEAプロット
library(enrichplot)
gseaplot2(gsea_result_comparison, geneSetID = 1)  # 上位のパスウェイを可視化


# ERKパスウェイの遺伝子を取得
erk_genes <- keggGet("hsa04010")[[1]]$GENE

# PI3K/Aktパスウェイの遺伝子を取得
pi3k_akt_genes <- keggGet("hsa04151")[[1]]$GENE

# MAPKシグナル伝達経路（hsa04010）から遺伝子データを取得
mapk_pathway <- keggGet("hsa04010")[[1]]$GENE

# 遺伝子データは、エントリが1つのリストとして返されるため、データを整形
mapk_genes <- as.character(mapk_pathway[seq(1, length(mapk_pathway), 2)])

# JNK関連の遺伝子（MAPK8, MAPK9, MAPK10など）をフィルタリング
jnk_genes <- grep("MAPK8|MAPK9|MAPK10", mapk_genes, value = TRUE)


# ERKシグナル伝達経路の遺伝子リストを拡充
erk_genes <- c("MAPK1", "MAPK3", "RAF1", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2")

# JNKシグナル伝達経路の遺伝子リスト
jnk_genes <- c("MAPK8", "MAPK9", "MAPK10", "JNK1", "JNK2", "JNK3")

# PI3K/Aktシグナル伝達経路の正しい遺伝子リスト
pi3k_akt_genes <- c("PIK3CA", "PIK3CB", "AKT1", "AKT2", "PTEN", "MTOR", "FOXO1", "FOXO3", "EGFR", "TGFA", "FGF1", "FGF2")

# ERK, JNK, PI3K/Akt の遺伝子セットをまとめる
phosphorylation_pathways <- data.frame(
  term = c(rep("ERK_PATHWAY", length(erk_genes)),
           rep("JNK_PATHWAY", length(jnk_genes)),
           rep("PI3K_AKT_PATHWAY", length(pi3k_akt_genes))),
  gene = c(erk_genes, jnk_genes, pi3k_akt_genes)
)

# GSEAの再実行
gsea_result <- GSEA(ranked_genes_comparison, TERM2GENE = phosphorylation_pathways, pvalueCutoff = 0.9, minGSSize = 1)

# 結果の確認と可視化
dotplot(gsea_result) + ggtitle("GSEA with ERK, JNK, and PI3K/Akt Pathways")

# mTORシグナル伝達経路の遺伝子リスト
mtor_genes <- c("MTOR", "RPTOR", "RICTOR", "AKT1", "AKT2", "PTEN")

# ERK, JNK, PI3K/Akt, mTOR の遺伝子セットをまとめる
phosphorylation_pathways <- data.frame(
  term = c(rep("ERK_PATHWAY", length(erk_genes)),
           rep("JNK_PATHWAY", length(jnk_genes)),
           rep("PI3K_AKT_PATHWAY", length(pi3k_akt_genes)),
           rep("MTOR_PATHWAY", length(mtor_genes))),
  gene = c(erk_genes, jnk_genes, pi3k_akt_genes, mtor_genes)
)

# GSEAの実行
gsea_result <- GSEA(ranked_genes_comparison, TERM2GENE = phosphorylation_pathways, pvalueCutoff = 0.9, minGSSize = 1)

# 結果の可視化
dotplot(gsea_result) + ggtitle("GSEA with ERK, JNK, PI3K/Akt, and mTOR Pathways")

# PI3K/AktとmTOR経路の遺伝子セットの確認
print(phosphorylation_pathways[phosphorylation_pathways$term == "PI3K_AKT_PATHWAY", ])
print(phosphorylation_pathways[phosphorylation_pathways$term == "MTOR_PATHWAY", ])
# PI3K/AktやmTOR経路に関連する遺伝子の発現量を確認
VlnPlot(object = data_IP_sampled, features = c("PIK3CA", "AKT1", "MTOR", "PTEN"))
# PI3K/AktやmTOR関連遺伝子のランキングリストでの位置を確認
ranked_genes_comparison[names(ranked_genes_comparison) %in% c("PIK3CA", "AKT1", "MTOR", "PTEN")]
# 特定のクラスターをサブセット化
subset_data <- subset(data_IP_sampled, idents = "Fibroblasts")

# 再度GSEAを実行
gsea_result <- GSEA(ranked_genes_comparison, TERM2GENE = phosphorylation_pathways, pvalueCutoff = 0.9, minGSSize = 1)

# 結果の確認と可視化
dotplot(gsea_result) + ggtitle("GSEA with ERK, JNK, PI3K/Akt, and mTOR Pathways (Subset)")



# ERKシグナル伝達経路の遺伝子リストを拡充
erk_genes <- c("MAPK1", "MAPK3", "RAF1", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2")

# JNKシグナル伝達経路の遺伝子リスト
jnk_genes <- c("MAPK8", "MAPK9", "MAPK10", "JNK1", "JNK2", "JNK3")

# PI3K/Aktシグナル伝達経路の遺伝子リスト
pi3k_akt_genes <- c("PIK3CA", "PIK3CB", "AKT1", "AKT2", "PTEN", "MTOR", "FOXO1", "FOXO3", "EGFR", "TGFA", "FGF1", "FGF2")

# mTORシグナル伝達経路の遺伝子リスト
mtor_genes <- c("MTOR", "RPTOR", "RICTOR", "AKT1", "AKT2", "PTEN")

# Wntシグナル伝達経路の遺伝子リスト
wnt_genes <- c("CTNNB1", "WNT1", "WNT3A", "FZD1", "LRP5", "LRP6")

# TGFβシグナル伝達経路の遺伝子リスト
tgfb_genes <- c("TGFB1", "TGFB2", "SMAD2", "SMAD3", "SMAD4", "SMAD7")

# NFκBシグナル伝達経路の遺伝子リスト
nfkb_genes <- c("NFKB1", "RELA", "IKBKG", "CHUK", "IKBKB")

# MAPK/MEKシグナル伝達経路の遺伝子リスト
mapk_mek_genes <- c("MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK14")

# Notchシグナル伝達経路の遺伝子リスト
notch_genes <- c("NOTCH1", "NOTCH2", "DLL1", "JAG1")

# 全シグナル伝達経路の遺伝子セットをまとめる
phosphorylation_pathways <- data.frame(
  term = c(rep("ERK_PATHWAY", length(erk_genes)),
           rep("JNK_PATHWAY", length(jnk_genes)),
           rep("PI3K_AKT_PATHWAY", length(pi3k_akt_genes)),
           rep("MTOR_PATHWAY", length(mtor_genes)),
           rep("WNT_PATHWAY", length(wnt_genes)),
           rep("TGF_BETA_PATHWAY", length(tgfb_genes)),
           rep("NFKB_PATHWAY", length(nfkb_genes)),
           rep("MAPK_MEK_PATHWAY", length(mapk_mek_genes)),
           rep("NOTCH_PATHWAY", length(notch_genes))),
  gene = c(erk_genes, jnk_genes, pi3k_akt_genes, mtor_genes, wnt_genes, tgfb_genes, nfkb_genes, mapk_mek_genes, notch_genes)
)

# GSEAの実行
gsea_result <- GSEA(ranked_genes_comparison, TERM2GENE = phosphorylation_pathways, pvalueCutoff = 0.9, minGSSize = 1)

# GSEA結果の確認
head(gsea_result@result)

# 結果の可視化
dotplot(gsea_result) + ggtitle("GSEA with Multiple Signaling Pathways")

#やり直し
fibroblasts_data$LRP1_group <- ifelse(fibroblasts_data@assays$SCT@data["LRP1", ] > median(fibroblasts_data@assays$SCT@data["LRP1", ]), "LRP1_high", "LRP1_low")
Idents(fibroblasts_data) <- "LRP1_group"
DEGs_LRP1 <- FindMarkers(fibroblasts_data, ident.1 = "LRP1_high", ident.2 = "LRP1_low", logfc.threshold = 0.000001)
ranked_genes_LRP1 <- DEGs_LRP1$avg_log2FC
names(ranked_genes_LRP1) <- rownames(DEGs_LRP1)
ranked_genes_LRP1 <- sort(ranked_genes_LRP1, decreasing = TRUE)
gsea_result_LRP1 <- GSEA(ranked_genes_LRP1, TERM2GENE = pathway_db, pvalueCutoff = 1, minGSSize = 10, maxGSSize = 1000)
gsea_plot <- gseaplot2(gsea_result_LRP1, geneSetID = "KEGG_MAPK_SIGNALING_PATHWAY")
print(gsea_plot)
# シグナル伝達経路の遺伝子リストを拡充
erk_genes <- c("MAPK1", "MAPK3", "RAF1", "BRAF", "MEK1", "MEK2", "ERK1", "ERK2")
jnk_genes <- c("MAPK8", "MAPK9", "MAPK10", "JNK1", "JNK2", "JNK3")
pi3k_akt_genes <- c("PIK3CA", "PIK3CB", "AKT1", "AKT2", "PTEN", "MTOR", "FOXO1", "FOXO3", "EGFR", "TGFA", "FGF1", "FGF2")
mtor_genes <- c("MTOR", "RPTOR", "RICTOR", "AKT1", "AKT2", "PTEN")
wnt_genes <- c("CTNNB1", "WNT1", "WNT3A", "FZD1", "LRP5", "LRP6")
tgfb_genes <- c("TGFB1", "TGFB2", "SMAD2", "SMAD3", "SMAD4", "SMAD7")
nfkb_genes <- c("NFKB1", "RELA", "IKBKG", "CHUK", "IKBKB")
mapk_mek_genes <- c("MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "MAPK14")
notch_genes <- c("NOTCH1", "NOTCH2", "DLL1", "JAG1")

# 全シグナル伝達経路の遺伝子セットをまとめる
phosphorylation_pathways <- data.frame(
  term = c(rep("ERK_PATHWAY", length(erk_genes)),
           rep("JNK_PATHWAY", length(jnk_genes)),
           rep("PI3K_AKT_PATHWAY", length(pi3k_akt_genes)),
           rep("MTOR_PATHWAY", length(mtor_genes)),
           rep("WNT_PATHWAY", length(wnt_genes)),
           rep("TGF_BETA_PATHWAY", length(tgfb_genes)),
           rep("NFKB_PATHWAY", length(nfkb_genes)),
           rep("MAPK_MEK_PATHWAY", length(mapk_mek_genes)),
           rep("NOTCH_PATHWAY", length(notch_genes))),
  gene = c(erk_genes, jnk_genes, pi3k_akt_genes, mtor_genes, wnt_genes, tgfb_genes, nfkb_genes, mapk_mek_genes, notch_genes)
)

# GSEAの実行
gsea_result <- GSEA(ranked_genes_comparison, TERM2GENE = phosphorylation_pathways, pvalueCutoff = 0.9, minGSSize = 1)
head(gsea_result@result)
dotplot(gsea_result) + ggtitle("GSEA with Multiple Signaling Pathways")


# LRP1_highとLRP1_lowにおける特定の代謝経路のUMAPプロット
# LRP1の発現に基づいて高発現群と低発現群を分ける
fibroblasts_metabolism$LRP1_group <- ifelse(fibroblasts_metabolism@assays$SCT@data["LRP1", ] > median(fibroblasts_metabolism@assays$SCT@data["LRP1", ]), "LRP1_high", "LRP1_low")
# DotPlotで複数のパスウェイを比較
input.pathway <- c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(fibroblasts_metabolism, pathway = input.pathway, phenotype = "LRP1_group", norm = "y")

p <- DimPlot.metabolism(fibroblasts_metabolism, pathway = "Glycolysis / Gluconeogenesis", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1.5)
print(p)
# DotPlotで複数のパスウェイを比較
input.pathway <- c("Glycolysis / Gluconeogenesis", "Oxidative phosphorylation", "Citrate cycle (TCA cycle)")
DotPlot.metabolism(fibroblasts_metabolism, pathway = input.pathway, phenotype = "LRP1_group", norm = "y")

# 追加するパスウェイ
input.pathway <- c(
  "Glycolysis / Gluconeogenesis", 
  "Oxidative phosphorylation", 
  "Citrate cycle (TCA cycle)", 
  "Pentose phosphate pathway",            # 追加のパスウェイ
  "Fatty acid metabolism",                # 追加のパスウェイ
  "Purine metabolism",                    # 追加のパスウェイ
  "Pyrimidine metabolism",                # 追加のパスウェイ
  "Steroid biosynthesis"                  # 追加のパスウェイ
)

# DotPlotでパスウェイを比較
DotPlot.metabolism(fibroblasts_metabolism, pathway = input.pathway, phenotype = "LRP1_group", norm = "y")
# DotPlotでパスウェイを比較し、縦横を逆にする
DotPlot.metabolism(fibroblasts_metabolism, pathway = input.pathway, phenotype = "LRP1_group", norm = "y") + coord_flip()

#2024年11月5日追加解析

# データ準備
celltype_counts <- table(data_BLM@meta.data$celltype)
celltype_df <- as.data.frame(celltype_counts) %>%
  mutate(Percentage = Freq / sum(Freq) * 100)

ggplot(celltype_df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Composition in data_BLM (細分化された表示)",
       x = NULL, y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_discrete(name = "Cell Type")

ggplot(celltype_df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Composition in data_BLM (細分化された表示)",
       x = NULL, y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(name = "Cell Type", palette = "Set3")

ggplot(celltype_df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Composition in data_BLM (細分化された表示)",
       x = NULL, y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_viridis_d(name = "Cell Type", option = "magma")

ggplot(celltype_df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Composition in data_BLM (グレースケール表示)",
       x = NULL, y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_grey(name = "Cell Type", start = 0, end = 0.9)  # グレースケール範囲を指定

custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(celltype_df$Var1)))

ggplot(celltype_df, aes(x = "", y = Percentage, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Composition in data_BLM (細分化された表示)",
       x = NULL, y = "Percentage (%)") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(name = "Cell Type", values = custom_colors)

#２０２４年１２月１１日日付毎の遺伝子発現の経過をグラフ化
#Fibroblastsのオブジェクト作成
fibroblasts <- subset(data_BLM, subset = celltype %in% c("Fibroblast", "Fibroblast_Peribronchiolar"))

#SCTransformの実行
fibroblasts <- SCTransform(fibroblasts, verbose = FALSE)

#conditionから日付情報だけを抽出
fibroblasts@meta.data$day <- gsub("(day\\d+).*", "\\1", fibroblasts@meta.data$condition)

#日付順に並べ替え
fibroblasts@meta.data$day <- factor(fibroblasts@meta.data$day, 
                                    levels = c("day00", "day03", "day07", "day14", "day28", "day42", "day63"))

head(rownames(fibroblasts))

#遺伝子リスト
genes <- c("Lrp1", "Crp", "Sftpd", "Col1a1")

#遺伝子発現量を日付ごとにまとめる
expr_list <- lapply(genes, function(gene) {
  expr <- FetchData(fibroblasts, vars = c(gene, "day"))
  colnames(expr) <- c("expression", "day")
  expr_summary <- expr %>%
    group_by(day) %>%
    summarise(
      mean_expression = mean(expression),
      se_expression = sd(expression) / sqrt(n())
    )
  expr_summary$gene <- gene
  return(expr_summary)
})

#データを結合
plot_data <- do.call(rbind, expr_list)

#折れ線グラフの作成
# 正規化: 遺伝子ごとに発現量をスケーリング
plot_data <- plot_data %>%
  group_by(gene) %>%
  mutate(scaled_expression = scale(mean_expression))

# day14 を除外
plot_data_filtered <- plot_data %>%
  dplyr::filter(day != "day14")  # day14 を除外

# スケーリングされたデータをプロット
ggplot(plot_data_filtered, aes(x = day, y = scaled_expression, group = gene, color = gene)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = scaled_expression - se_expression, ymax = scaled_expression + se_expression), width = 0.2) +
  theme_minimal() +
  labs(
    title = "Gene Expression Trends Over Time (Scaled)",
    x = "Time Points",
    y = "Scaled Gene Expression",
    color = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#２０２４年１２月１１日ssGSEA解析を行い時系列毎のパスウェイ変化を図示
# 1. パスウェイデータの準備（指定されたGMTファイルを使用）
gmt_file <- "C:/Users/yamas/OneDrive/デスクトップ/RStudioの資料そのほか/c2.all.v2024.1.Hs.symbols.gmt"
gene_sets <- getGmt(gmt_file)

# 2. 遺伝子発現データの抽出
# fibroblastsオブジェクトから発現マトリックスを取得
expr_matrix <- GetAssayData(fibroblasts, slot = "data")

# 3. ssGSEA解析の実行
ssgsea_scores <- gsva(as.matrix(expr_matrix), gene_sets, method = "ssgsea", mx.diff = TRUE, verbose = TRUE)

# 4. パスウェイスコアを日付情報でまとめる
# スコアをデータフレーム化
scores_df <- as.data.frame(t(ssgsea_scores))  # パスウェイ×サンプル
scores_df$day <- fibroblasts@meta.data$day   # 日付情報を追加

# 平均と標準誤差を計算
scores_summary <- scores_df %>%
  group_by(day) %>%
  summarise(across(where(is.numeric), list(mean = ~ mean(.), se = ~ sd(.) / sqrt(n()))))

# 5. 可視化するデータの整形
# 選択したパスウェイを含むデータをピボットテーブル化
plot_data <- scores_summary %>% pivot_longer(cols = -day, names_to = c("pathway", ".value"), names_pattern = "(.+)_(mean|se)")

# 可視化するパスウェイを選択（例: 上位3つのパスウェイ）
selected_pathways <- c("CHUANG_OXIDATIVE_STRESS_RESPONSE_UP", "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON", "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_BLUE_UP")

# 選択したパスウェイのみフィルタリング
plot_data_filtered <- plot_data %>%
  dplyr::filter(pathway %in% selected_pathways) %>%
  group_by(day, pathway) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    se_score = sd(score, na.rm = TRUE) / sqrt(n())
  )


# 6. ggplotで時系列グラフを作成
ggplot(plot_data_filtered, aes(x = day, y = mean_score, group = pathway, color = pathway)) +
  geom_line(size = 1) +  # 折れ線グラフ
  geom_point(size = 3) +  # データポイント
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score), width = 0.2) +  # エラーバー
  theme_minimal() +
  labs(
    title = "Pathway Activity Trends Over Time (ssGSEA)",
    x = "Time Points",
    y = "Pathway Activity Score (mean ± SE)",
    color = "Pathway"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # x軸ラベルを傾ける
    legend.title = element_text(size = 10),  # 凡例タイトルサイズ
    legend.text = element_text(size = 9)    # 凡例テキストサイズ
  )


# 3. パスウェイスコアの統計解析
# データを長い形式に変換
ssgsea_long <- scores_df %>% pivot_longer(cols = -day, names_to = "pathway", values_to = "score")

# パスウェイごとにANOVAを実行して時系列の有意差を評価
anova_results <- ssgsea_long %>%
  group_by(pathway) %>%
  summarise(p_value = summary(aov(score ~ day, data = cur_data()))[[1]][["Pr(>F)"]][1]) %>%
  arrange(p_value)  # p値でソート

# 有意なパスウェイ（例: p < 0.05）を抽出
significant_pathways <- anova_results %>% dplyr::filter(p_value < 0.05) %>% pull(pathway)

# 亢進パスウェイの確認
print("Significant Pathways:")
print(significant_pathways)

# 4. 有意なパスウェイを可視化
# 有意なパスウェイのスコアをフィルタリング
plot_data <- ssgsea_long %>% dplyr::filter(pathway %in% significant_pathways)  # dplyr::filter を明示的に指定


# プロット
ggplot(plot_data, aes(x = day, y = score, group = pathway, color = pathway)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    title = "Significantly Enriched Pathways Over Time (ssGSEA)",
    x = "Time Points",
    y = "Pathway Activity Score",
    color = "Pathway"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#２０２５年７月１０日 時系列毎にLRP1の発現推移グラフ　PROGENY解析

# 対象細胞種を定義
selected_celltypes <- c(
  "Fibro_proliferated", "Fibro_peribron", "Fibro_alveolar", "Fibro_adven",
  "Mac_AM", "Mac_IM", "Mac_int", "Mac_proliferated"
)

# 対象細胞を抽出
subset_data <- subset(data_BLM, 
                      subset = celltype_coarse %in% selected_celltypes)

# Lrp1の発現値を抽出
subset_data$Lrp1_expr <- FetchData(subset_data, vars = "Lrp1")[, 1]

# 共通化された DayGroup を作成
subset_data$DayGroup <- str_extract(subset_data$condition, "^day\\d+")

# 平均発現を CellType × DayGroup ごとに集計（df_avgを再定義！）
df_avg <- subset_data@meta.data %>%
  group_by(celltype_coarse, DayGroup) %>%
  summarise(mean_expr = mean(Lrp1_expr), .groups = "drop") %>%
  dplyr::rename(CellType = celltype_coarse, Day = DayGroup)

# day56除外 → Dayをfactorに設定 → Zスコア変換（この順！）
df_avg_filtered <- df_avg %>%
  dplyr::filter(Day %in% c("day00", "day03", "day07", "day14", "day28", "day42")) %>%
  mutate(Day = factor(Day, levels = c("day00", "day03", "day07", "day14", "day28", "day42"))) %>%
  group_by(CellType) %>%
  mutate(zscore_expr = scale(mean_expr)[, 1]) %>%
  ungroup()

# プロット
ggplot(df_avg_filtered, aes(x = Day, y = CellType, fill = zscore_expr)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Z-score\nLrp1") +
  theme_minimal(base_size = 14) +
  labs(title = "Lrp1 Expression (Z-score)", x = "Day", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#PROGENY解析
# 対象細胞種を定義
selected_celltypes <- c(
  "Fibro_proliferated", "Fibro_peribron", "Fibro_alveolar", "Fibro_adven",
  "Mac_AM", "Mac_IM", "Mac_int", "Mac_proliferated"
)

# 抽出
subset_data <- subset(data_BLM, subset = celltype_coarse %in% selected_celltypes)

# 平均発現行列（log-normalized）
avg_expr <- AverageExpression(
  subset_data,
  group.by = c("celltype_coarse", "condition"),
  return.seurat = FALSE
)$RNA


progeny_scores <- progeny(
  expr = as.matrix(avg_expr),
  scale = TRUE,        # log-normalized データでOK
  organism = "Mouse",
  top = 500
)

# 使用できるパスウェイ名を確認
head(colnames(progeny_scores))
# → 例: "Androgen" "EGFR" "Estrogen" "Hypoxia" "JAK-STAT" "MAPK"

df_progeny <- melt(t(progeny_scores))
colnames(df_progeny) <- c("Pathway", "Group", "Score")

# CellTypeとDayに分解（"CellType_condition" から）
df_progeny <- df_progeny %>%
  dplyr::mutate(
    Day = str_extract(Group, "day\\d+"),
    CellType = sub("_day\\d+.*", "", Group)
  )

# 不要なday除去・順序設定
df_progeny <- df_progeny %>%
  dplyr::filter(Day %in% c("day00", "day03", "day07", "day14", "day28", "day42")) %>%
  dplyr::mutate(Day = factor(Day, levels = c("day00", "day03", "day07", "day14", "day28", "day42")))

ggplot(dplyr::filter(df_progeny, Pathway == "TGFb"),
       aes(x = Day, y = CellType, fill = Score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = median(df_progeny$Score, na.rm = TRUE),
                       name = "PROGENy Score") +
  theme_minimal(base_size = 14) +
  labs(title = "TGFβ Pathway Activity", x = "Day", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# 例：複数パスウェイをfacetで並べて表示
selected_pathways <- c("TGFb", "MAPK", "NFkB", "EGFR", "Hypoxia", "JAK-STAT", "PI3K", "TGFb", "TNFa")

ggplot(df_progeny %>% dplyr::filter(Pathway %in% selected_pathways),
       aes(x = Day, y = CellType, fill = Score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, name = "PROGENy Score") +
  facet_wrap(~Pathway) +
  theme_minimal(base_size = 14) +
  labs(title = "Pathway Activities", x = "Day", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#CellchatとNicheNet解析の追加
# 対象細胞タイプ（Mac_AMが送信、Fibro_系が受信）
celltypes_of_interest <- c("Mac_AM", "Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar")
# 対象日付（必要に応じて増減可能）
data_BLM$condition_merged <- case_when(
  grepl("day00", data_BLM$condition) ~ "day00",
  grepl("day03", data_BLM$condition) ~ "day03",
  grepl("day07", data_BLM$condition) ~ "day07",
  grepl("day14", data_BLM$condition) ~ "day14",
  grepl("day28", data_BLM$condition) ~ "day28",
  grepl("day42", data_BLM$condition) ~ "day42",
  grepl("day63", data_BLM$condition) ~ "day63",
  TRUE ~ "other"
)

day_points <- c("day00", "day03", "day07", "day14", "day28")

DefaultAssay(data_BLM)
SeuratObject::Assays(data_BLM)
data_BLM[["RNA"]]@counts[1:5, 1:5]

data.input <- GetAssayData(obj_day, assay = "RNA", slot = "counts")
meta.data <- obj_day@meta.data

cellchat <- createCellChat(object = data.input, meta = meta.data, group.by = "celltype_coarse")

cellchat_list <- list()

for (day in day_points) {
  cat("Processing", day, "\n")
  
  obj_day <- subset(data_BLM, subset = condition_merged == day & celltype_coarse %in% celltypes_of_interest)
  
  data.input <- GetAssayData(obj_day, assay = "RNA", slot = "counts")
  meta.input <- obj_day@meta.data
  
  cellchat <- createCellChat(object = data.input, meta = meta.input, group.by = "celltype_coarse")
  cellchat@DB <- CellChatDB.mouse
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  
  cellchat <- filterCommunication(cellchat)  # ← この時点で warning が出るのはOK（Fibro_proliferatedが除外されただけ）
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  # 保存
  cellchat_list[[day]] <- cellchat
}  # ←



# day00の例（TGFβ signaling）
netVisual_chord_gene(
  object = cellchat_list[["day00"]],
  signaling = "TGFb",
  sources.use = "Mac_AM",
  targets.use = c("Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar"),
  lab.cex = 0.9
)

# 貢献度確認（netAnalysis_contribution）
netAnalysis_contribution(cellchat_list[["day00"]], signaling = "TGFb")

plot_list <- lapply(day_points, function(day) {
  netVisual_chord_gene(cellchat_list[[day]],
                       signaling = "TGFb",
                       sources.use = "Mac_AM",
                       targets.use = c("Fibro_proliferated", "Fibro_peribron", "Fibro_adven", "Fibro_alveolar"),
                       title.name = paste("TGFb -", day))
})

wrap_plots(plotlist = plot_list, ncol = 2)

netVisual_aggregate(cellchat_list[[day]], signaling = "TGFb", layout = "circle")

netVisual_chord_gene(cellchat_list[[day]],
                     signaling = "TGFb",
                     sources.use = NULL,  # すべての送信元を対象
                     targets.use = NULL,  # すべての受信先を対象
                     title.name = paste("TGFb -", day))


#PROGENYとLRP1発現相関
# 行：CellType、列：Day の wide型データフレームが前提
cor_results <- expand.grid(CellType = unique(df_avg_filtered$CellType),
                           Pathway = unique(df_progeny$Pathway))


cor_results$correlation <- mapply(function(ct, pw){
  lrp1 <- df_avg_filtered %>% dplyr::filter(CellType == ct) %>% arrange(Day) %>% pull(zscore_expr)
  prog <- df_progeny %>% dplyr::filter(CellType == ct, Pathway == pw) %>% arrange(Day) %>% pull(Score)
  
  # デバッグ出力
  cat("CellType:", ct, " Pathway:", pw, 
      " LRP1 Length:", length(lrp1), 
      " PROGENy Length:", length(prog), "\n")
  
  if (length(lrp1) == length(prog) & length(lrp1) > 1) {
    stats::cor(lrp1, prog, method = "pearson")
  } else {
    NA
  }
}, cor_results$CellType, cor_results$Pathway)

unique(df_avg_filtered$Day)
unique(df_progeny$Day)

df_progeny_unique <- df_progeny %>% 
  dplyr::distinct(CellType, Pathway, Day, .keep_all = TRUE)

cor_results$correlation <- mapply(function(ct, pw){
  lrp1 <- df_avg_filtered %>% dplyr::filter(CellType == ct) %>% arrange(Day) %>% pull(zscore_expr)
  prog <- df_progeny_unique %>% dplyr::filter(CellType == ct, Pathway == pw) %>% arrange(Day) %>% pull(Score)
  
  if (length(lrp1) == length(prog) & length(lrp1) > 1) {
    stats::cor(lrp1, prog, method = "pearson")
  } else {
    NA
  }
}, cor_results$CellType, cor_results$Pathway)

ggplot(cor_results, aes(x = Pathway, y = CellType, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                       name = "Correlation") +
  theme_minimal(base_size = 14) +
  labs(title = "Correlation between LRP1 and PROGENy Pathways",
       x = "Pathway", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Nichenet
# Seurat IDを celltype_coarse に設定
Idents(data_BLM) <- "celltype_coarse"

# 1. 送信側: Mac_AM（day14のみに絞る）
sender_all <- subset(data_BLM, subset = celltype_coarse == "Mac_AM" & condition_merged == "day14")
sender_all$Lrp1_expr <- FetchData(sender_all, vars = "Lrp1")[, 1]
highLRP1_sender <- subset(sender_all, subset = Lrp1_expr > quantile(Lrp1_expr, 0.75))

# 2. 受信側: Fibro群（day14のみに絞る）
receiver_cells <- subset(data_BLM, subset = celltype_coarse %in% c("Fibro_proliferated", "Fibro_alveolar", "Fibro_adven", "Fibro_peribron") & condition_merged == "day14")

# 3. 発現行列の取得
exprs_sender <- GetAssayData(highLRP1_sender, slot = "data")  # log-normalized
exprs_receiver <- GetAssayData(receiver_cells, slot = "data")

# 4. ターゲット遺伝子セットの定義（線維化マーカー）
geneset_oi <- c("Col1a1", "Col3a1", "Fn1", "Acta2", "Tgfb1")

# 送信側（Mac_AM, day14 の LRP1 high）から発現している遺伝子を抽出
sender_genes <- rownames(exprs_sender)[Matrix::rowSums(exprs_sender > 0.1) / ncol(exprs_sender) > 0.10]

# 受信側から発現している遺伝子を抽出（背景遺伝子として使う）
receiver_genes <- rownames(exprs_receiver)[Matrix::rowSums(exprs_receiver > 0.1) / ncol(exprs_receiver) > 0.10]

# NicheNetの既知リガンド一覧（ligand_target_matrixとlr_networkから取得）
ligands_in_network <- lr_network %>% pull(from) %>% unique()

# 送信側で発現していて、かつネットワークに含まれるリガンドの交差を取る
potential_ligands <- intersect(sender_genes, ligands_in_network)

# 必要なデータの読み込み（初回だけでOK）
organism <- "mouse"
ligand_target_matrix <- readRDS("C:/Users/yamas/OneDrive/デスクトップ/RStudioの資料そのほか/ligand_target_matrix_nsga2r_final_mouse.rds")

# ligand-receptor network # weighted integrated ligand–target matrix
lr_network <- readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks <- readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# SeuratのFindMarkersなどで DEGs を抽出した後
deg_result <- FindMarkers(data_BLM,
                          ident.1 = "day14",
                          ident.2 = "day00",
                          group.by = "condition_merged",
                          subset.ident = "Fibro_proliferated")

# Ensemblのマウス遺伝子名変換
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# gene名 → Ensemblまたは別の形式へ
converted <- getBM(
  filters = "mgi_symbol",
  attributes = c("mgi_symbol", "external_gene_name", "ensembl_gene_id", "gene_biotype"),
  values = deg_result$gene,
  mart = mart
)

geneset_oi <- converted %>%
  filter(external_gene_name %in% rownames(ligand_target_matrix)) %>%
  pull(external_gene_name) %>%
  unique()


intersect(geneset_oi, rownames(ligand_target_matrix))

# リガンドアクティビティスコアの計算
ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = receiver_genes,
  potential_ligands = potential_ligands,
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network
)


# 上位リガンドを表示
head(ligand_activities[order(-ligand_activities$aupr), ], 10)

# 可視化（棒グラフ）
top_ligands <- ligand_activities %>% top_n(10, wt = aupr)
ggplot(top_ligands, aes(x = reorder(test_ligand, aupr), y = aupr)) +
  geom_col() +
  coord_flip() +
  labs(x = "Ligand", y = "Activity (AUPR)", title = "Top predicted ligands affecting fibroblast fibrosis genes") +
  theme_minimal()

VlnPlot(data_BLM, features = c("F7", "C3"), group.by = "celltype", split.by = "condition_merged")

# F7, C3 に対応するターゲット遺伝子（上位のもの）を抽出
ligands_of_interest <- c("F7", "C3")

# ligand_target_matrix の中で非ゼロのターゲット遺伝子を抽出
target_prediction_df <- expand.grid(
  ligand = ligands_of_interest,
  target = rownames(ligand_target_matrix)
)

# 各リガンド・ターゲット組み合わせに対してスコア（weight）を取得
target_prediction_df$weight <- mapply(function(lig, tar) {
  ligand_target_matrix[tar, lig]
}, target_prediction_df$ligand, target_prediction_df$target)

# 重み順に上位のターゲット遺伝子を表示
top_targets <- target_prediction_df %>%
  filter(weight > 0) %>%
  arrange(desc(weight)) %>%
  group_by(ligand) %>%
  slice_head(n = 20)

# 表示
print(top_targets)


# リガンド・ターゲット予測スコア取得（すでに行っていればスキップ）
nichenet_get_weighted_target_predictions <- function(ligand_target_matrix, ligands){
  ligand_target_matrix = ligand_target_matrix[intersect(rownames(ligand_target_matrix), rownames(ligand_target_matrix)), , drop = FALSE]
  ligand_target_df = ligand_target_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("target") %>% 
    pivot_longer(-target, names_to = "ligand", values_to = "weight") %>% 
    filter(ligand %in% ligands) %>%
    filter(!is.na(weight)) %>%
    arrange(-weight)
  return(ligand_target_df)
}


ligand_target_df <- nichenet_get_weighted_target_predictions(
  ligand_target_matrix = ligand_target_matrix,
  ligands = c("F7", "C3")
)

# 表示したいリガンド・ターゲットの上位（必要に応じて絞る）
top_targets <- ligand_target_df %>%
  group_by(ligand) %>%
  top_n(5, wt = weight) %>%  # 上位5ターゲットなど
  ungroup()

# circos plotの作成
nichenet_circle_ligand_target_plot <- function(
    ligand_target_df,
    geneset,
    background_genes,
    title = "Ligand-target circos plot",
    cutoff = 0.01,
    color_ligand = "red",
    color_target = "steelblue",
    color_background = "grey90"
){
  library(circlize)
  library(dplyr)
  
  circos.clear()
  
  ligand_target_df <- ligand_target_df %>% 
    filter(target %in% geneset | target %in% background_genes)
  
  ligands <- unique(ligand_target_df$ligand)
  targets <- unique(ligand_target_df$target)
  
  sectors <- c(ligands, targets)
  factors <- factor(sectors, levels = sectors)
  circos.par(start.degree = 90, gap.degree = 3)
  circos.initialize(factors = factors, xlim = cbind(rep(0, length(factors)), rep(1, length(factors))))
  
  circos.trackPlotRegion(
    ylim = c(0, 1), panel.fun = function(x, y) {
      sector_name = get.cell.meta.data("sector.index")
      circos.text(CELL_META$xcenter, 0.5, sector_name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
    },
    bg.border = NA
  )
  
  for(i in 1:nrow(ligand_target_df)) {
    ligand <- ligand_target_df$ligand[i]
    target <- ligand_target_df$target[i]
    weight <- ligand_target_df$weight[i]
    
    if (weight > cutoff) {
      col <- ifelse(target %in% geneset, color_target, color_background)
      circos.link(ligand, runif(1), target, runif(1), col = col, border = NA)
    }
  }
  
  title(title, cex.main = 1.2)
}


nichenet_circle_ligand_target_plot(
  ligand_target_df = top_targets,
  geneset = geneset_oi,  # 自分が注目している fibrosis genesなど
  background_genes = receiver_genes,
  title = "Ligand-target circos plot"
)


#Connectome解析
day_list <- SplitObject(data_BLM, split.by = "condition_merged")

connectome_list <- list()

for (day in names(day_list)) {
  seurat_day <- day_list[[day]]
  
  # 必要に応じて NormalizeData と ScaleData（SCTransformedであれば不要）
  seurat_day <- NormalizeData(seurat_day)
  seurat_day <- ScaleData(seurat_day)
  
  # Connectome の作成
  connectome_list[[day]] <- CreateConnectome(
    object = seurat_day,
    species = "mouse",
    group.by = "celltype",  # ここは cell type 情報のカラム名
    min.cells = 5,
    p.values = TRUE  # P値も計算（時間がかかる）
  )
}

# 使用する細胞種の指定
# 例: day00のサブセットを作成して Connectome オブジェクトに変換
data_day00 <- subset(data_BLM, subset = condition_merged == "day00")

connectome_day00 <- CreateConnectome(
  data_day00,
  species = "mouse",  # humanの場合は "human"
  p.values = TRUE
)
cells.of.interest <- c('Mac_AM','Fibro_proliferated','Fibro_peribron','Fibro_adven','Fibro_alveolar')

CircosPlot(connectome_day00,
           sources.use = cells.of.interest,
           targets.use = cells.of.interest,
           min.pct = 0.1,
           min.z = 1,
           cols.use = NULL,
           lab.cex = 0.6,
           title = "Connectome Circos Plot - day00")

for (day in names(connectome_list)) {
  message("Plotting: ", day)
  CircosPlot(connectome_list[[day]],
             sources.use = cells.of.interest,
             targets.use = cells.of.interest,
             min.pct = 0.1,
             min.z = 1,
             lab.cex = 0.6,
             title = paste("Connectome Circos Plot -", day))
}

CircosPlot(connectome_day00,
           sources.include = cells.of.interest,
           targets.include = cells.of.interest,
           min.pct = 0.1,
           min.z = 1,
           cols.use = NULL,
           lab.cex = 0.6,
           title = "Connectome Circos Plot - day00")



cells.of.interest <- c('Mac_AM','Fibro_proliferated','Fibro_peribron','Fibro_Adven','Fibro_alveolar')
CircosPlot(connectome_day00,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from normalized slot')
CircosPlot(data_BLM,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,lab.cex = 0.6,title = 'Edgeweights from scaled slot')
CircosPlot(data_BLM,weight.attribute = 'weight_norm',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from normalized slot)')
CircosPlot(data_BLM,weight.attribute = 'weight_sc',sources.include = cells.of.interest,targets.include = cells.of.interest,balanced.edges = F,lab.cex = 0.6,title = 'Ligand vs. receptor expression (from scaled slot)')
CircosPlot(data_BLM,sources.include = 'Mac_AM',lab.cex = 0.6)
