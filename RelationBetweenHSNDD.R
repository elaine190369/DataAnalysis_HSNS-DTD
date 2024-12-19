print('start')
if (!requireNamespace("CCA", quietly = TRUE)) install.packages("CCA")
if (!requireNamespace("qgraph", quietly = TRUE)) install.packages("qgraph")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("psych", quietly = TRUE)) install.packages("psych")
if (!requireNamespace("NetworkComparisonTest", quietly = TRUE)) install.packages("NetworkComparisonTest")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

library(CCA)
library(qgraph)
library(igraph)
library(psych)
library(NetworkComparisonTest)
library(Matrix)

# 读取
raw_data <- read.table("data.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
raw_data <- na.omit(raw_data)
raw_data <- raw_data[-1, ]
raw_data <- as.data.frame(lapply(raw_data[, 1:22], as.numeric))
print(dim(raw_data))
str(raw_data)
summary(raw_data)

raw_data <- raw_data[!rowSums(raw_data == 0) > 0, ]
data1 <- raw_data[, 1:10]
data2 <- raw_data[, 11:22]

# 检查数据
print(dim(data1)) 
print(dim(data2)) 
str(data1)
str(data2)
summary(data1)
summary(data2)

#区分效度：皮尔逊相关
raw_data$totalHSN <- rowSums(raw_data[, 1:10])
raw_data$totalDTD <- rowSums(raw_data[, 11:22])
print('this is pearson')
print(cor.test(raw_data$totalHSN, raw_data$totalDTD, method = "pearson"))
data <- raw_data[, 1:22]

# 1. 执行典型相关分析
cca_result <- cancor(data1, data2)

# 输出结果
print("典型相关系数：")
print(cca_result$cor)

# 变量的典型权重（载荷系数）
print("量表1的典型变量载荷系数：")
print(cca_result$xcoef)
print("量表2的典型变量载荷系数：")
print(cca_result$ycoef)
print("变量集合 X 的权重（第一个典型变量 U1）：")
print(cca_result$xcoef[, 1])
print("变量集合 Y 的权重（第一个典型变量 V1）：")
print(cca_result$ycoef[, 1])

# 2. 建立网络
# 计算相关矩阵
cor_matrix <- cor(data, use = "pairwise.complete.obs")  # 计算相关矩阵
cor_matrix <- as.matrix(nearPD(cor_matrix)$mat)        # 逼近为正定矩阵
glasso <- EBICglasso(cor_matrix, n = nrow(raw_data), gamma = 0.10)
glasso <- abs(glasso)

# 绘制网络图
png("output_plots_Rela/NetworkBetweenHSNandDTD.png", width = 800, height = 600)
qgraph(glasso, layout = "spring", labels = colnames(data), groups = list(
  "Scale1" = 1:ncol(data1),
  "Scale2" = (ncol(data1) + 1):(ncol(data1) + ncol(data2))
))
dev.off()

# qgraph(glasso, layout = "spring", title = "Full Sample Network", theme = "colorblind", labels = colnames(raw_data), vsize = 10, esize = 10)
# dev.off()

print('startBridge')
# 桥接强度
bridge_results <- bridge(glasso)
summary(bridge_results)

# 可视化桥接中心性指标
png("output_plots_Rela/Bridge_Centrality.png", width = 800, height = 600)
print(plot(bridge_results, include = c("Bridge Strength", "Bridge Betweenness", "Bridge Closeness")))
dev.off()
print('network finished')

# 6. 网络稳定性评估
print(dim(data)) 
str(data)
summary(data)

boot_results <- bootnet(
  data,                   
  type = "person",
  nBoots = 1000,
  default = 'EBICglasso',
  statistics = c("edge", "strength", "betweenness", "closeness", "distance")
)

# 检查稳定性结果
png("output_plots_Rela/BootResults.png", width = 800, height = 600)
print(plot(boot_results))
dev.off()
print('done')
