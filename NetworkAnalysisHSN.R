if (!require("qgraph")) install.packages("qgraph", dependencies = TRUE) 
if (!require("bootnet")) install.packages("bootnet", dependencies = TRUE) 
if (!require("NetworkComparisonTest")) install.packages("NetworkComparisonTest", dependencies = TRUE) 
if (!require("networktools")) install.packages("networktools", dependencies = TRUE) 
if (!require("reshape2")) install.packages("reshape2", dependencies = TRUE) 

# 加载必要的库
library(qgraph)
library(bootnet)
library(NetworkComparisonTest)
library(networktools)
library(reshape2)

# 1. 读取
raw_data <- read.table("data.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data <- raw_data[, 1:10]
data <- na.omit(data)
data <- data[-1, ]
data <- as.data.frame(lapply(data, as.numeric))
data <- data[!rowSums(data == 0) > 0, ]
# 检查数据
print(dim(data)) 
str(data)
summary(data)


# 计算皮尔逊相关矩阵
cor_matrix <- cor(data, use = "pairwise.complete.obs")
print(cor_matrix)

# 2. 分析 HSN 的网络结构
# 使用 EBICglasso 方法估计网络
glasso_full <- EBICglasso(cor_matrix, n = nrow(data), gamma = 0.10)

# 保存
png("output_plots_HSNnetwork/Network.png", width = 800, height = 800)
qgraph(glasso_full, layout = "spring", title = "Full Sample Network", theme = "colorblind", labels = colnames(data[, 1:10]), vsize = 10, esize = 10)
dev.off()

# matrix_graph <- as.matrix(network$graph)
# # 检查网络边权重是否有负值
# summary(matrix_graph)
# # 如果存在负边权重，将其设置为绝对值
# matrix_graph[matrix_graph < 0] <- abs(matrix_graph)
# print(matrix_graph)

# 3. 网络比较检验
#分性别数据
raw_data <- as.data.frame(lapply(raw_data, as.numeric))
print(class(raw_data))
print(colnames(raw_data))
data2 <- raw_data[ ,c(1:10, 24)]
data2 <- na.omit(data2) 
data2 <- data2[!rowSums(data == 0) > 0, ]
group1_data <- data2[data2$gender == 1, ]
group2_data <- data2[data2$gender == 2, ]
str(group1_data)
str(group2_data)


# 计算两个组的相关矩阵，并将负边界专为绝对值
cor_matrix_group1 <- cor(group1_data[, 1:10], use = "pairwise.complete.obs")
cor_matrix_group2 <- cor(group2_data[, 1:10], use = "pairwise.complete.obs")
print(cor_matrix_group1)
print(cor_matrix_group2)

#负权重不能参与，通过检测替换为绝对值
cor_matrix_group1 <- abs(cor_matrix_group1)
cor_matrix_group2 <- abs(cor_matrix_group2)

# 使用 EBICglasso 方法估计网络
glasso1 <- EBICglasso(cor_matrix_group1, n = nrow(group1_data[, 1:10]), gamma = 0.10)
png("output_plots_HSNnetwork/MaleNetwork.png", width = 800, height = 800)
qgraph(glasso_full, layout = "spring", title = "Male Sample Network", theme = "colorblind", labels = colnames(group1_data[, 1:10]), vsize = 10, esize = 10)
dev.off()

glasso2 <- EBICglasso(cor_matrix_group2, n = nrow(group2_data[, 1:10]), gamma = 0.10)
png("output_plots_HSNnetwork/FemaleNetwork.png", width = 800, height = 800)
qgraph(glasso_full, layout = "spring", title = "Female Sample Network", theme = "colorblind", labels = colnames(group2_data[, 1:10]), vsize = 10, esize = 10)
dev.off()

# 运行NCT
print("startNCT!!!")
nct_test <- NCT(data1 = group1_data[, 1:10], 
                data2 = group2_data[, 1:10],
                gamma = 0.50, 
                test.edges = TRUE, 
                edges = "all",
                abs=TRUE,
                test.centrality = TRUE, 
                p.adjust.methods = "holm",
                paired = FALSE, 
                it = 1000)

summary(nct_test)

print('this is graph1')
png("output_plots_HSNnetwork/ComparisonEdgeDiff.png", width = 800, height = 800)
# 计算边权重差异
edge_diff <- glasso1 - glasso2

color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks <- seq(min(edge_diff, na.rm = TRUE), max(edge_diff, na.rm = TRUE), length.out = 101)

# 绘制热图
image(1:nrow(edge_diff), 1:ncol(edge_diff), t(edge_diff[nrow(edge_diff):1, ]), 
      col = color_palette, breaks = breaks, xlab = "Node 1", ylab = "Node 2", 
      main = "Edge Weight Difference", axes = FALSE)

# 添加轴标签
axis(1, at = 1:nrow(edge_diff), labels = rownames(edge_diff), las = 2, cex.axis = 0.7)
axis(2, at = 1:ncol(edge_diff), labels = colnames(edge_diff), las = 2, cex.axis = 0.7)

# 添加图例
legend("topright", legend = round(seq(min(edge_diff, na.rm = TRUE), 
                                      max(edge_diff, na.rm = TRUE), length.out = 5), 2),
       fill = color_palette[c(1, 25, 50, 75, 100)], title = "Edge\nDifference",
       cex = 0.8, bty = "n")

# 关闭图形设备
dev.off()


# 4. 计算中心性指标
# 使用 bootnet 包计算中心性指标
glasso_full <- abs(glasso_full)
centrality <- centralityPlot(glasso_full, include = c("Strength", "Closeness", "Betweenness"))

# 提取中心性指标
centrality_indices <- centrality_auto(glasso_full)
print(centrality_indices)
print('ok')

# 5. 桥接中心性分析
# 使用 networktools 包计算桥接中心性
bridge_results <- bridge(glasso_full)
summary(bridge_results)
print('yes')

# 可视化桥接中心性指标
png("output_plots_HSNnetwork/Bridge_Centrality.png", width = 800, height = 600)
print(plot(bridge_results, include = c("Bridge Strength", "Bridge Betweenness", "Bridge Closeness")))
dev.off()
print('good')

# 提取中心性指标
# 使用 bootnet 包分别计算 glasso1 和 glasso2 的中心性指标
glasso1 <- abs(glasso1)
glasso2 <- abs(glasso2)

centrality1 <- centralityPlot(glasso1, include = c("Strength", "Closeness", "Betweenness"), verbose = FALSE)
centrality2 <- centralityPlot(glasso2, include = c("Strength", "Closeness", "Betweenness"), verbose = FALSE)

# 提取中心性指标数据
centrality_indices1 <- centrality_auto(glasso1)  # glasso1 的中心性指标
centrality_indices2 <- centrality_auto(glasso2)  # glasso2 的中心性指标

# 打印中心性指标
print('is it ok?')
print(centrality_indices1)
print(centrality_indices2)
print('photo now')

# 提取桥接中心性指标
# 提取桥接中心性指标
strength1 <- centrality_indices1$node.centrality[, "Strength"]
betweenness1 <- centrality_indices1$node.centrality[, "Betweenness"]
closeness1 <- centrality_indices1$node.centrality[, "Closeness"]

strength2 <- centrality_indices2$node.centrality[, "Strength"]
betweenness2 <- centrality_indices2$node.centrality[, "Betweenness"]
closeness2 <- centrality_indices2$node.centrality[, "Closeness"]

# 节点名称
nodes <- rownames(centrality_indices1$node.centrality)  # 从中心性指标矩阵中获取节点名称
if (is.null(nodes)) {
  nodes <- paste("Node", 1:nrow(glasso1))  # 如果没有名称，则创建默认名称
}

# 检查数据长度
print("Lengths check:")
print(paste("Nodes length:", length(nodes)))
print(paste("Strength1 length:", length(strength1)))
print(paste("Strength2 length:", length(strength2)))
print(paste("Betweenness1 length:", length(betweenness1)))
print(paste("Betweenness2 length:", length(betweenness2)))
print(paste("Closeness1 length:", length(closeness1)))
print(paste("Closeness2 length:", length(closeness2)))

# 检查是否有NA值
print("NA check:")
print(sum(is.na(strength1)))
print(sum(is.na(strength2)))
print(sum(is.na(betweenness1)))
print(sum(is.na(betweenness2)))
print(sum(is.na(closeness1)))
print(sum(is.na(closeness2)))

# 绘图函数
plot_data <- function(x1, x2, nodes, title, xlab) {
  # 确保所有向量长度一致
  n <- length(nodes)
  if(length(x1) != n || length(x2) != n) {
    warning(paste("Length mismatch in", title))
    return(FALSE)
  }
  
  # 移除NA值
  valid <- !is.na(x1) & !is.na(x2)
  if(!all(valid)) {
    warning(paste("NA values found in", title))
    x1 <- x1[valid]
    x2 <- x2[valid]
    nodes <- nodes[valid]
  }
  
  # 绘图
  plot(
    x1, 1:length(nodes), type = "o", col = "blue", pch = 16, cex = 1.2,
    xlim = range(c(x1, x2), na.rm = TRUE),
    ylim = c(1, length(nodes)),
    xlab = xlab, ylab = "", main = title, axes = FALSE
  )
  lines(x2, 1:length(nodes), type = "o", col = "red", pch = 16, cex = 1.2)
  axis(1)
  axis(2, at = 1:length(nodes), labels = nodes, las = 2, cex.axis = 0.8)
  box()
  legend("topright", legend = c("glasso1", "glasso2"), 
         col = c("blue", "red"), lty = 1, pch = 16, bty = "n", cex = 1)
  
  return(TRUE)
}

# 绘制图形
png("output_plots_HSNnetwork/Centrality_Comparison.png", width = 1800, height = 800)
par(mfrow = c(1, 3), mar = c(5, 8, 4, 2))

# 绘制三个图
plot_data(strength1, strength2, nodes, "Strength", "Strength")
plot_data(betweenness1, betweenness2, nodes, "Betweenness", "Betweenness")
plot_data(closeness1, closeness2, nodes, "Closeness", "Closeness")

dev.off()

# 6. 网络稳定性评估
# 使用 bootnet 包评估网络模型稳定性
boot_results <- bootnet(data, 
                        type = "person",
                        nBoots = 1000, # 引导抽样次数
                        default = "EBICglasso",
                        statistics = c("edge", "strength", "betweenness", "closeness", "distance")) # 使用 EBICglasso

# 检查稳定性结果
png("output_plots_HSNnetwork/BootResults.png", width = 800, height = 600)
print(plot(boot_results))
dev.off()
print('done')

