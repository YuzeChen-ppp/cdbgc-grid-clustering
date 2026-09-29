# ==============================================================================
# 📊 多算法对决：DBCV网格密度 vs K-means(K=2) vs 层次聚类(K=2) (ARI 对比版)
# ==============================================================================
required_packages <- c("ggplot2", "igraph", "patchwork", "dbscan")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(ggplot2)
library(igraph)
library(patchwork)
library(dbscan)

# ==========================================
# 0. 数据集构建：生成带有随机噪声的双月牙形
# ==========================================
set.seed(42) 
n_samples <- 200

cluster1 <- data.frame(
  x = cos(theta1 <- seq(0, pi, length.out = n_samples)) * 3 + rnorm(n_samples, 0, 0.25),
  y = sin(theta1) * 3 + rnorm(n_samples, 0, 0.25),
  true_label = "Original Cluster 1"
)

cluster2 <- data.frame(
  x = 3 - cos(theta2 <- seq(0, pi, length.out = n_samples)) * 3 + rnorm(n_samples, 0, 0.25),
  y = 1.5 - sin(theta2) * 3 + rnorm(n_samples, 0, 0.25),
  true_label = "Original Cluster 2"
)

outliers <- data.frame(
  x = runif(15, -4, 7),  
  y = runif(15, -4, 5),  
  true_label = "True Outlier"
)

data <- rbind(cluster1, cluster2, outliers)
true_vector <- as.numeric(factor(data$true_label)) # 真实标签向量

# ==========================================
# 1. 算法 I：带 Outlier 惩罚的网格密度聚类 (自适应识别)
# ==========================================
run_clustering <- function(input_data, grid_n, thr) {
  x_breaks <- seq(min(input_data$x), max(input_data$x), length.out = grid_n + 1)
  y_breaks <- seq(min(input_data$y), max(input_data$y), length.out = grid_n + 1)
  ix <- findInterval(input_data$x, x_breaks, all.inside = TRUE)
  iy <- findInterval(input_data$y, y_breaks, all.inside = TRUE)
  Lx <- x_breaks[ix]; Rx <- x_breaks[ix + 1]
  Ly <- y_breaks[iy]; Ry <- y_breaks[iy + 1]
  tx <- (input_data$x - Lx) / (Rx - Lx); ty <- (input_data$y - Ly) / (Ry - Ly)
  Wx <- cbind((1 - tx)/2, 0.5, tx/2); Wy <- cbind((1 - ty)/2, 0.5, ty/2)
  
  density_vec <- numeric(grid_n * grid_n)
  for (a in 1:3) {
    for (b in 1:3) {
      xi <- ix + (a - 2); yj <- iy + (b - 2); w  <- Wx[, a] * Wy[, b]
      valid <- xi >= 1 & xi <= grid_n & yj >= 1 & yj <= grid_n
      if (any(valid)) {
        flat_idx <- xi[valid] + (yj[valid] - 1) * grid_n
        chunk_sum <- rowsum(w[valid], flat_idx)
        density_vec[as.integer(rownames(chunk_sum))] <- density_vec[as.integer(rownames(chunk_sum))] + chunk_sum[, 1]
      }
    }
  }
  density_mat <- matrix(density_vec, nrow = grid_n, ncol = grid_n)
  active <- density_mat >= thr; coords <- which(active, arr.ind = TRUE)
  
  cluster_grid <- matrix(NA, grid_n, grid_n)
  if (nrow(coords) > 0) {
    adj_matrix <- as.matrix(dist(coords, method = "maximum")) <= 1
    diag(adj_matrix) <- 0
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    cluster_grid[coords] <- components(g)$membership
  }
  
  unique_clusters <- sort(unique(cluster_grid[!is.na(cluster_grid)]))
  if (length(unique_clusters) > 0) {
    point_cluster_weight <- matrix(0, nrow = nrow(input_data), ncol = length(unique_clusters))
    for (a in 1:3) {
      for (b in 1:3) {
        xi <- ix + (a - 2); yj <- iy + (b - 2); valid <- xi >= 1 & xi <= grid_n & yj >= 1 & yj <= grid_n
        c_id <- rep(NA, nrow(input_data))
        if (any(valid)) c_id[valid] <- cluster_grid[cbind(xi[valid], yj[valid])]
        w <- Wx[, a] * Wy[, b]
        for (k_idx in seq_along(unique_clusters)) {
          k <- unique_clusters[k_idx]; matches_k <- !is.na(c_id) & c_id == k
          point_cluster_weight[matches_k, k_idx] = point_cluster_weight[matches_k, k_idx] + w[matches_k]
        }
      }
    }
    max_idx <- apply(point_cluster_weight, 1, which.max)
    res_cluster <- unique_clusters[max_idx]
    res_cluster[point_cluster_weight[cbind(1:nrow(input_data), max_idx)] <= 1/4] <- -1 
  } else { res_cluster <- rep(-1, nrow(input_data)) }
  return(list(cluster = res_cluster))
}

# 网格搜索最优参数
grid_n_settings <- seq(10, 40, by = 2); thr_settings <- seq(0.2, 1.5, by = 0.1)   
tuning_grid <- expand.grid(grid_n = grid_n_settings, thr = thr_settings)
tuning_grid$DBCV <- NA; outlier_penalty_weight <- 2.5 

for (i in 1:nrow(tuning_grid)) {
  res <- run_clustering(data, tuning_grid$grid_n[i], tuning_grid$thr[i])
  cl_labels <- res$cluster; cl_labels[cl_labels == -1] <- 0
  num_cl <- length(unique(cl_labels[cl_labels != 0]))
  if (num_cl < 2 || num_cl > 5 || any(table(cl_labels[cl_labels != 0]) < 5)) { tuning_grid$DBCV[i] <- -2.0; next }
  raw_dbcv <- tryCatch({
    dbcv_out <- dbscan::dbcv(as.matrix(data[, c("x", "y")]), cl_labels)
    if (is.list(dbcv_out)) ifelse(!is.null(dbcv_out$dbcv), dbcv_out$dbcv, dbcv_out[[1]]) else dbcv_out
  }, error = function(e) { -1.0 })
  tuning_grid$DBCV[i] <- raw_dbcv - (outlier_penalty_weight * (sum(cl_labels == 0) / length(cl_labels)))
}
best_idx <- which.max(tuning_grid$DBCV)
best_res <- run_clustering(data, tuning_grid$grid_n[best_idx], tuning_grid$thr[best_idx])

# 计算网格密度的 ARI
data$cluster_density <- best_res$cluster
ari_density <- igraph::compare(true_vector, as.numeric(factor(data$cluster_density)), method = "adjusted.rand")

# ==========================================
# 2. 🎯 算法 II：K-means 聚类 (强制设定 K = 2)
# ==========================================
km_fit <- kmeans(data[, c("x", "y")], centers = 2, nstart = 25)
data$cluster_kmeans <- km_fit$cluster

# 计算 K-means 的 ARI
ari_kmeans <- igraph::compare(true_vector, data$cluster_kmeans, method = "adjusted.rand")

# ==========================================
# 3. 🎯 算法 III：层次聚类 (强制设定 K = 2)
# ==========================================
distance_matrix <- dist(data[, c("x", "y")])
hc_fit <- hclust(distance_matrix, method = "ward.D2") 
data$cluster_hc <- cutree(hc_fit, k = 2)

# 计算层次聚类的 ARI
ari_hc <- igraph::compare(true_vector, data$cluster_hc, method = "adjusted.rand")

# ==========================================
# 4. 📈 输出指标对比报告
# ==========================================
report_df <- data.frame(
  Algorithm = c("Grid Density", "K-means (Forced K=2)", "Hierarchical (Forced K=2)"),
  Target_Clusters = c("Auto Detect (2 + Outliers)", "Strictly 2", "Strictly 2"),
  ARI_Score = c(ari_density, ari_kmeans, ari_hc)
)

cat("\n========================================================\n")
cat("📊  算法外部指标评估报告 (Adjusted Rand Index) \n")
cat("========================================================\n")
print(report_df, row.names = FALSE)
cat("========================================================\n\n")

# ==========================================
# 5. 🎨 终极全要素对比看板渲染
# ==========================================
theme_set(theme_minimal() + theme(panel.grid = element_blank(), plot.title = element_text(face="bold", size=12)))

# 统一配色方案 (1和2表示聚类，-1或0表示网格识别出的噪点)
palette_colors <- c("1"="#377eb8", "2"="#4daf4a", "-1"="black")

p_density <- ggplot(data, aes(x, y, color = factor(cluster_density))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = palette_colors) +
  labs(title = "1. Grid Density", subtitle = paste("ARI:", round(ari_density, 4)), color = "Cluster")

p_kmeans <- ggplot(data, aes(x, y, color = factor(cluster_kmeans))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = palette_colors) +
  labs(title = "2. K-means (K=2)", subtitle = paste("ARI:", round(ari_kmeans, 4)), color = "Cluster")

p_hc <- ggplot(data, aes(x, y, color = factor(cluster_hc))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = palette_colors) +
  labs(title = "3. Hierarchical (K=2)", subtitle = paste("ARI:", round(ari_hc, 4)), color = "Cluster")

# 使用 patchwork 渲染分栏看板
final_comparison_board <- p_density + p_kmeans + p_hc + plot_layout(ncol = 3, guides = "collect")
print(final_comparison_board)                    