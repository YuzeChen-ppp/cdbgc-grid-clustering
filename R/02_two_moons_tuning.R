
required_packages <- c("ggplot2", "igraph", "patchwork", "dbscan")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(ggplot2)
library(igraph)
library(patchwork)
library(dbscan)

# ==========================================
# 0. 数据集构建：生成带有随机噪声的双月牙形（Two Moons）
# ==========================================
set.seed(42) 
n_samples <- 200

# 月牙 1
theta1 <- seq(0, pi, length.out = n_samples)
cluster1 <- data.frame(
  x = cos(theta1) * 3 + rnorm(n_samples, 0, 0.25),
  y = sin(theta1) * 3 + rnorm(n_samples, 0, 0.25),
  true_label = "Original Cluster 1"
)

# 月牙 2
theta2 <- seq(0, pi, length.out = n_samples)
cluster2 <- data.frame(
  x = 3 - cos(theta2) * 3 + rnorm(n_samples, 0, 0.25),
  y = 1.5 - sin(theta2) * 3 + rnorm(n_samples, 0, 0.25),
  true_label = "Original Cluster 2"
)

# 背景随机噪声（Outliers）
outliers <- data.frame(
  x = runif(15, -4, 7),  
  y = runif(15, -4, 5),  
  true_label = "True Outlier"
)

# 合并完整数据集
data <- rbind(cluster1, cluster2, outliers)

# ==========================================
# 1. 核心算法功能函数
# ==========================================
run_clustering <- function(input_data, grid_n, thr) {
  x_breaks <- seq(min(input_data$x), max(input_data$x), length.out = grid_n + 1)
  y_breaks <- seq(min(input_data$y), max(input_data$y), length.out = grid_n + 1)
  
  ix <- findInterval(input_data$x, x_breaks, all.inside = TRUE)
  iy <- findInterval(input_data$y, y_breaks, all.inside = TRUE)
  
  Lx <- x_breaks[ix]; Rx <- x_breaks[ix + 1]
  Ly <- y_breaks[iy]; Ry <- y_breaks[iy + 1]
  tx <- (input_data$x - Lx) / (Rx - Lx)
  ty <- (input_data$y - Ly) / (Ry - Ly)
  
  Wx <- cbind((1 - tx)/2, 0.5, tx/2)
  Wy <- cbind((1 - ty)/2, 0.5, ty/2)
  
  density_vec <- numeric(grid_n * grid_n)
  for (a in 1:3) {
    for (b in 1:3) {
      xi <- ix + (a - 2)
      yj <- iy + (b - 2)
      w  <- Wx[, a] * Wy[, b]
      valid <- xi >= 1 & xi <= grid_n & yj >= 1 & yj <= grid_n
      if (any(valid)) {
        flat_idx <- xi[valid] + (yj[valid] - 1) * grid_n
        chunk_sum <- rowsum(w[valid], flat_idx)
        idx_int   <- as.integer(rownames(chunk_sum))
        density_vec[idx_int] <- density_vec[idx_int] + chunk_sum[, 1]
      }
    }
  }
  density_mat <- matrix(density_vec, nrow = grid_n, ncol = grid_n)
  
  active <- density_mat >= thr
  coords <- which(active, arr.ind = TRUE)
  
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
        xi <- ix + (a - 2)
        yj <- iy + (b - 2)
        valid <- xi >= 1 & xi <= grid_n & yj >= 1 & yj <= grid_n
        c_id <- rep(NA, nrow(input_data))
        if (any(valid)) c_id[valid] <- cluster_grid[cbind(xi[valid], yj[valid])]
        w <- Wx[, a] * Wy[, b]
        for (k_idx in seq_along(unique_clusters)) {
          k <- unique_clusters[k_idx]
          matches_k <- !is.na(c_id) & c_id == k
          point_cluster_weight[matches_k, k_idx] = point_cluster_weight[matches_k, k_idx] + w[matches_k]
        }
      }
    }
    max_idx <- apply(point_cluster_weight, 1, which.max)
    max_w   <- point_cluster_weight[cbind(1:nrow(input_data), max_idx)]
    res_cluster <- unique_clusters[max_idx]
    res_cluster[max_w <= 1/4] <- -1 
  } else {
    res_cluster <- rep(-1, nrow(input_data))
  }
  
  return(list(cluster = res_cluster, cluster_grid = cluster_grid, density_vec = density_vec))
}

# ==========================================
# 2.Outlier 惩罚权重的 DBCV 
# ==========================================
grid_n_settings <- seq(10, 40, by = 2)       
thr_settings    <- seq(0.2, 1.5, by = 0.1)   

tuning_grid <- expand.grid(grid_n = grid_n_settings, thr = thr_settings)
tuning_grid$DBCV <- NA
tuning_grid$num_clusters <- NA


outlier_penalty_weight <- 2.5 

for (i in 1:nrow(tuning_grid)) {
  res <- run_clustering(data, tuning_grid$grid_n[i], tuning_grid$thr[i])
  cl_labels <- res$cluster
  
  #  1. DBCV 规范对齐
  cl_labels[cl_labels == -1] <- 0
  
  #  2. 基础拓扑合法性前置检查（防止 C++ 崩溃）
  valid_points_mask <- cl_labels != 0
  valid_clusters    <- unique(cl_labels[valid_points_mask])
  num_cl            <- length(valid_clusters)
  tuning_grid$num_clusters[i] <- num_cl
  
  if (num_cl < 2 || num_cl > 5) {
    tuning_grid$DBCV[i] <- -2.0 # 连基础大骨架都不符合，直接判定极低分
    next
  }
  cluster_sizes <- table(cl_labels[valid_points_mask])
  if (any(cluster_sizes < 5)) {
    tuning_grid$DBCV[i] <- -2.0 # 包含极极端碎屑簇
    next
  }
  
  #  3. 安全计算标准 DBCV 分数
  raw_dbcv <- tryCatch({
    dbcv_out <- dbscan::dbcv(as.matrix(data[, c("x", "y")]), cl_labels)
    if (is.list(dbcv_out)) {
      if (!is.null(dbcv_out$dbcv)) dbcv_out$dbcv else dbcv_out[[1]]
    } else {
      dbcv_out
    }
  }, error = function(e) {
    -1.0 
  })
  

  # 计算当前配置下，全场被当成噪声（0）的点所占的比例
  noise_ratio <- sum(cl_labels == 0) / length(cl_labels)
  
  # 最终得分 = 拓扑密集度得分 - 噪声比例惩罚项
  # 如果强行把点变成 Outlier 来让原始 DBCV 变高，减法惩罚项会直接让最终指标发生崩塌
  tuning_grid$DBCV[i] <- raw_dbcv - (outlier_penalty_weight * noise_ratio)
}

#  Outlier 惩罚约束后的全球最优解 
best_idx    <- which.max(tuning_grid$DBCV)
best_grid_n <- tuning_grid$grid_n[best_idx]
best_thr    <- tuning_grid$thr[best_idx]


# ==========================================
# 3. 提取最优解实体指标
# ==========================================
best_res     <- run_clustering(data, best_grid_n, best_thr)
data$cluster <- best_res$cluster
cluster_grid <- best_res$cluster_grid
density_vec  <- best_res$density_vec

x_breaks <- seq(min(data$x), max(data$x), length.out = best_grid_n + 1)
y_breaks <- seq(min(data$y), max(data$y), length.out = best_grid_n + 1)

# ==========================================
# 4. 渲染
# ==========================================
cluster_ids <- sort(unique(data$cluster[data$cluster != -1]))
cols        <- setNames(rainbow(length(cluster_ids)), cluster_ids)
true_cols   <- c("Original Cluster 1" = "#377eb8", "Original Cluster 2" = "#4daf4a", "True Outlier" = "#e41a1c")

# ---  图 1 DBCV 2D 超参数平面 ---
p_tune_2d <- ggplot(tuning_grid, aes(x = grid_n, y = thr)) +
  geom_tile(aes(fill = DBCV), color = "gray95", lwd = 0.05) +
  scale_fill_viridis_c(option = "viridis") +
  geom_point(data = tuning_grid[best_idx, ], aes(x = grid_n, y = thr), 
             color = "red", size = 4, shape = 18) + 
  theme_minimal() +
  labs(
    title = "Plot 1: Outlier-Penalized DBCV Space",
    subtitle = paste("Best spot forced by penalty: grid_n =", best_grid_n, ", thr =", best_thr),
    x = "Grid Resolution (grid_n)", y = "Density Threshold (thr)", fill = "Score"
  )

# ---  图 2：空间密度分布 ---
grid_centers <- expand.grid(ix = 1:best_grid_n, iy = 1:best_grid_n)
grid_centers$x <- (x_breaks[grid_centers$ix] + x_breaks[grid_centers$ix + 1]) / 2
grid_centers$y <- (y_breaks[grid_centers$iy] + y_breaks[grid_centers$iy + 1]) / 2
grid_centers$density <- density_vec

p2 <- ggplot(grid_centers, aes(x, y)) +
  geom_tile(aes(fill = density), color = "white", lwd = 0.2) +
  scale_fill_viridis_c(option = "magma") + theme_minimal() +
  labs(fill = "Density", title = "Plot 2: Density Map") +
  theme(panel.grid = element_blank())

# ---🧩 图 3：过滤后的连通网格骨架 ---
grid_centers$grid_cluster <- cluster_grid[cbind(grid_centers$ix, grid_centers$iy)]

p3 <- ggplot(grid_centers, aes(x, y)) +
  geom_tile(aes(fill = factor(grid_cluster)), color = "white", lwd = 0.4) +
  scale_fill_manual(values = cols, na.value = "gray95") + theme_minimal() +
  labs(fill = "Grid Cluster", title = "Plot 3: Grid Subgraph Components") +
  theme(panel.grid = element_blank())

# ---  图 4 & 5：带最优方格线的分栏盲测结果对比 ---
p4 <- ggplot(data, aes(x, y)) +
  geom_vline(xintercept = x_breaks, color = "gray90", linetype = "dashed", lwd = 0.3) + 
  geom_hline(yintercept = y_breaks, color = "gray90", linetype = "dashed", lwd = 0.3) + 
  geom_point(aes(color = true_label), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = true_cols) + theme_minimal() + 
  labs(color = "Ground Truth", title = "Plot 4: Ground Truth", subtitle = paste("Grid Overlay:", best_grid_n, "x", best_grid_n)) +
  theme(panel.grid = element_blank())

p1 <- ggplot(data, aes(x, y)) +
  geom_vline(xintercept = x_breaks, color = "gray90", linetype = "dashed", lwd = 0.3) + 
  geom_hline(yintercept = y_breaks, color = "gray90", linetype = "dashed", lwd = 0.3) + 
  geom_point(aes(color = factor(cluster)), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(cols, "-1" = "black")) + theme_minimal() +
  labs(color = "Calculated", title = "Plot 5: Heavily Penalized Output", subtitle = paste("Outliers (Black) minimized due to high penalty weight")) +
  theme(panel.grid = element_blank())

comparison_plot <- p4 + p1 + plot_layout(ncol = 2)


print(p_tune_2d)       
print(p2)              
print(p3)              
print(comparison_plot)