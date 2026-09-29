# ==========================================
# 0. Data (加入原始身份标签，含 (0,0) 噪点)
# ==========================================
set.seed(123)

cluster1 <- data.frame(x = rnorm(45, 0, 0.6), y = rnorm(45, 0, 0.6))
cluster1$true_label <- "Original Cluster 1"

cluster2 <- data.frame(x = rnorm(45, 5, 0.6), y = rnorm(45, 5, 0.6))
cluster2$true_label <- "Original Cluster 2"

outliers <- data.frame(x = runif(10, -3, 8),  y = runif(10, -3, 8))
outliers <- rbind(outliers, data.frame(x = 0, y = 0))
outliers$true_label <- "True Outlier"

data     <- rbind(cluster1, cluster2, outliers)

# ==========================================
# 1. Grid (10x10)
# ==========================================
grid_n <- 10
x_breaks <- seq(min(data$x), max(data$x), length.out = grid_n + 1)
y_breaks <- seq(min(data$y), max(data$y), length.out = grid_n + 1)

# ==========================================
# 2. 向量化权重计算与累加
# ==========================================
ix <- findInterval(data$x, x_breaks, all.inside = TRUE)
iy <- findInterval(data$y, y_breaks, all.inside = TRUE)

Lx <- x_breaks[ix]; Rx <- x_breaks[ix + 1]
Ly <- y_breaks[iy]; Ry <- y_breaks[iy + 1]

tx <- (data$x - Lx) / (Rx - Lx)
ty <- (data$y - Ly) / (Ry - Ly)

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

# ==========================================
# 3. 向量化建图与连通域
# ==========================================
thr <- 1
active <- density_mat >= thr
coords <- which(active, arr.ind = TRUE)
n_nodes <- nrow(coords)

cluster_grid <- matrix(NA, grid_n, grid_n)

if (n_nodes > 0) {
  library(igraph)
  adj_matrix <- as.matrix(dist(coords, method = "maximum")) <= 1
  diag(adj_matrix) <- 0
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  comp <- components(g)
  cluster_grid[coords] <- comp$membership
}

# ==========================================
# 4. 软权重分配点所属聚类 (1/3 贡献率规则)
# ==========================================
unique_clusters <- sort(unique(cluster_grid[!is.na(cluster_grid)]))

if (length(unique_clusters) > 0) {
  point_cluster_weight <- matrix(0, nrow = nrow(data), ncol = length(unique_clusters))
  colnames(point_cluster_weight) <- unique_clusters
  
  for (a in 1:3) {
    for (b in 1:3) {
      xi <- ix + (a - 2)
      yj <- iy + (b - 2)
      valid <- xi >= 1 & xi <= grid_n & yj >= 1 & yj <= grid_n
      
      c_id <- rep(NA, nrow(data))
      if (any(valid)) {
        c_id[valid] <- cluster_grid[cbind(xi[valid], yj[valid])]
      }
      
      w <- Wx[, a] * Wy[, b]
      
      for (k_idx in seq_along(unique_clusters)) {
        k <- unique_clusters[k_idx]
        matches_k <- !is.na(c_id) & c_id == k
        point_cluster_weight[matches_k, k_idx] <- point_cluster_weight[matches_k, k_idx] + w[matches_k]
      }
    }
  }
  
  max_idx <- apply(point_cluster_weight, 1, which.max)
  max_w   <- point_cluster_weight[cbind(1:nrow(data), max_idx)]
  
  data$cluster <- unique_clusters[max_idx]
  data$cluster[max_w <= 1/3] <- -1
} else {
  data$cluster <- -1
}

# ==========================================
# 5. 独立输出与合并对比输出
# ==========================================
library(ggplot2)
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

# 基础颜色配置
cluster_ids <- sort(unique(data$cluster[data$cluster != -1]))
cols        <- setNames(rainbow(length(cluster_ids)), cluster_ids)
true_cols   <- c("Original Cluster 1" = "#377eb8", "Original Cluster 2" = "#4daf4a", "True Outlier" = "#e41a1c")

# --- 【Plot 4：原始真实标签（左图准备）】 ---
p4 <- ggplot(data, aes(x, y)) +
  geom_vline(xintercept = x_breaks, color = "gray85", linetype = "dashed") +
  geom_hline(yintercept = y_breaks, color = "gray85", linetype = "dashed") +
  geom_point(aes(color = true_label), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = true_cols) +
  theme_minimal() +
  labs(
    color = "Ground Truth", 
    title = "A. Original Data Labels",
    subtitle = "Red points are pure mathematical noise"
  ) +
  theme(panel.grid = element_blank())

# --- 【Plot 1：算法计算结果（右图准备）】 ---
p1 <- ggplot(data, aes(x, y)) +
  geom_vline(xintercept = x_breaks, color = "gray85", linetype = "dashed") +
  geom_hline(yintercept = y_breaks, color = "gray85", linetype = "dashed") +
  geom_point(aes(color = factor(cluster)), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c(cols, "-1" = "black")) +
  theme_minimal() +
  labs(
    color = "Calculated Cluster", 
    title = "B. Algorithm's Result",
    subtitle = "Black points are classified as noise"
  ) +
  theme(panel.grid = element_blank())

# --- 【中间辅助图：网格密度热力图】 ---
grid_centers <- expand.grid(ix = 1:grid_n, iy = 1:grid_n)
grid_centers$x <- (x_breaks[grid_centers$ix] + x_breaks[grid_centers$ix + 1]) / 2
grid_centers$y <- (y_breaks[grid_centers$iy] + y_breaks[grid_centers$iy + 1]) / 2
grid_centers$density <- density_vec

p2 <- ggplot(grid_centers, aes(x, y)) +
  geom_tile(aes(fill = density), color = "white", lwd = 0.3) +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(fill = "Density", title = "Plot 2: Grid Density Heatmap") +
  theme(panel.grid = element_blank())

# --- 【中间辅助图：网格连通块图】 ---
grid_centers$grid_cluster <- cluster_grid[cbind(grid_centers$ix, grid_centers$iy)]

p3 <- ggplot(grid_centers, aes(x, y)) +
  geom_tile(aes(fill = factor(grid_cluster)), color = "white", lwd = 0.5) +
  scale_fill_manual(values = cols, na.value = "gray95") +
  theme_minimal() +
  labs(fill = "Grid Cluster", title = "Plot 3: Grid Connected Components") +
  theme(panel.grid = element_blank())


# --- 📌 顺次打印输出 📌 ---
print(p2)  # 先看热力图
print(p3)  # 再看连通块图

# 最后输出核心对比图：将原始真实标签(p4)和算法结果(p1)拼在同一张图片里
comparison_plot <- p4 + p1 + plot_layout(ncol = 2)
print(comparison_plot)