# 加载所需 R 包
print('ready')
if (!require("psych")) install.packages("psych", dependencies = TRUE)  # 心理学分析工具包
if (!require("GPArotation")) install.packages("GPArotation", dependencies = TRUE) # 因子旋转
if (!require("factoextra")) install.packages("factoextra", dependencies = TRUE)  # 因子分析工具包
if (!require("ltm")) install.packages("ltm", dependencies = TRUE)  # 项目分析
if (!require("ggplot2")) install.packages("ggplot2") # 可视化
if (!require("semPlot")) install.packages("semPlot")

library(psych)
library(GPArotation)
library(lavaan)
library(ltm)
library(factoextra)
library(ggplot2)
library(semPlot)

# 1. 数据导入
# 读取数据
data <- read.table("HSNS+DD/data.csv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
data <- na.omit(data)
data <- data[, 11:22]
data <- data[-1, ]
data <- as.data.frame(lapply(data, as.numeric))
data <- data[!rowSums(data == 0) > 0, ]
# 检查数据
print(dim(data)) 
str(data)
summary(data)

# 2. 项目分析
# 2.1 描述每个项目的选项分布
item_distributions <- apply(data, 2, function(column){
  freq_table <- table(column)
  freq_table/length(column)
})
print("各项目的选项分布：")
print(item_distributions)

# 2.2 计算区分度（项目-总分相关系数）
data$Total <- rowSums(data)

# 极端分组法计算区分度
# 划分高分组和低分组
n <- nrow(data)
high_group <- data[order(-data$Total), ][1:ceiling(n * 0.27), ]
low_group <- data[order(data$Total), ][1:ceiling(n * 0.27), ]

# 计算每题的高低分组平均分差
distinction <- sapply(1:(ncol(data) - 1), function(i) {
  high_mean <- mean(high_group[[i]])
  low_mean <- mean(low_group[[i]])
  high_mean - low_mean
})

# 输出极端分组法的区分度
names(distinction) <- colnames(data)[1:(ncol(data) - 1)]
print("区分度:")
print(distinction)

# 相关法计算区分度
# 计算每个项目的项目-总分相关
total_corr <- sapply(1:(ncol(data) - 1), function(i){
  corrected_total <- data$Total 
  cor(data[[i]], corrected_total)           
})

item_total_corr <- sapply(1:(ncol(data) - 1), function(i){
  corrected_total <- data$Total - data[[i]]  # 删去当前项目后的总分
  cor(data[[i]], corrected_total)          
})

correlation_matrix <- sapply(1:(ncol(data) - 1), function(i) {
  corrected_total <- data$Total - data[[i]]
  sapply(1:(ncol(data) - 1), function(j) {
    cor(data[[j]], corrected_total)
  })
})

colnames(correlation_matrix) <- colnames(data)[1:(ncol(data) - 1)]
rownames(correlation_matrix) <- colnames(data)[1:(ncol(data) - 1)]

print("每个项目与未校正的总分的相关矩阵：")
print(total_corr)
print("每个项目与校正后的总分的相关矩阵：")
print(item_total_corr)
print("每个项目与校正对应项目后的总分的相关矩阵：")
print(correlation_matrix)

# 3. 信度分析 
# 3.1 内部一致性信度
cronbach_alpha <- psych::alpha(data)
print("Cronbach's α：")
print(cronbach_alpha$total[1])  # 提取 Cronbach's α 值

# 3.2 分半信度
split_half <- psych::splitHalf(data)
print("分半信度：")
print(split_half)

library(lavaan)

# 定义模型
model <- '
  ML1 =~ V11 + V12 + V13 + V14
  ML2 =~ V15 + V16 + V17 + V18 
  ML3 =~ V19 + V20 + V21 + V22
'

# 执行验证性因子分析
fit <- cfa(model, data = data)

# 查看模型结果
summary(fit, fit.measures = TRUE, standardized = TRUE)

# 绘制最终模型路径图
semPaths(fit, "std", layout = "tree", residuals = FALSE, whatLabels = "std")

# 输出最终 CFA 模型的完整结果
cat("\n--- Final CFA Model Summary ---\n")
summary(fit, fit.measures = TRUE, standardized = TRUE)

# 提取并打印最终模型的标准化载荷
cat("\n--- Final Standardized Factor Loadings ---\n")
final_loadings <- standardizedSolution(fit)
final_loadings <- final_loadings[final_loadings$op == "=~", c("lhs", "rhs", "est.std")]
print(final_loadings)

if (!require("qgraph")) install.packages("qgraph", dependencies = TRUE) 
if (!require("bootnet")) install.packages("bootnet", dependencies = TRUE) 
if (!require("NetworkComparisonTest")) install.packages("NetworkComparisonTest", dependencies = TRUE) 
if (!require("networktools")) install.packages("networktools", dependencies = TRUE) 
if (!require("tidyr")) install.packages("tidyr", dependencies = TRUE) 

# 加载必要的库
library(qgraph)
library(bootnet)
library(NetworkComparisonTest)
library(networktools)

# 1. 读取
raw_data <- read.table("HSNS+DD/data.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
data <- raw_data[, 11:22]
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

# 2. 分析 DTD 的网络结构
# 使用 EBICglasso 方法估计网络
glasso_full <- EBICglasso(cor_matrix, n = nrow(data), gamma = 0.10)

# 加载 qgraph 包
library(qgraph)

# 使用 EBICglasso 方法估计网络
network <- qgraph::EBICglasso(cor_matrix, n = nrow(data))

# 绘制网络
qgraph(network, layout = "spring", title = "Network Visualization")

# 使用 qgraph 提取中心性指标
centrality <- centrality(network)
print(centrality)

# 可视化中心性指标
centralityPlot(network)

qgraph(
  network,
  layout = "spring",
  nodeNames = colnames(data),  # 可选：节点名称
  color = "lightblue",         # 节点颜色
  title = "Estimated Network",
  vsize = centrality$Strength  # 节点大小基于强度中心性
)

png("output_plots_DTDnetwork/network_visualization.png", width = 800, height = 800)
qgraph(network, layout = "spring", title = "Network Visualization")
dev.off()

# 3. 网络比较检验
#分性别数据
raw_data <- as.data.frame(lapply(raw_data, as.numeric))
print(class(raw_data))
print(colnames(raw_data))
data2 <- raw_data[ ,c(11:22, 24)]
data2 <- na.omit(data2) 
data2 <- data2[!rowSums(data == 0) > 0, ]
group1_data <- data2[data2$gender == 1, ]
group2_data <- data2[data2$gender == 2, ]
str(group1_data)
str(group2_data)

# 计算两个组的相关矩阵，并将负边界专为绝对值
cor_matrix_group1 <- cor(group1_data[, 1:12], use = "pairwise.complete.obs")
cor_matrix_group2 <- cor(group2_data[, 1:12], use = "pairwise.complete.obs")
print(cor_matrix_group1)
print(cor_matrix_group2)

#负权重不能参与，通过检测替换为绝对值
cor_matrix_group1 <- abs(cor_matrix_group1)
cor_matrix_group2 <- abs(cor_matrix_group2)

# 使用 EBICglasso 方法估计网络
glasso1 <- EBICglasso(cor_matrix_group1, n = nrow(group1_data[, 1:12]), gamma = 0.10)
png("output_plots_DTDnetwork/MaleNetwork.png", width = 800, height = 800)
qgraph(glasso_full, layout = "spring", title = "Male Sample Network", theme = "colorblind", labels = colnames(group1_data[, 1:12]), vsize = 10, esize = 10)
dev.off()

glasso2 <- EBICglasso(cor_matrix_group2, n = nrow(group2_data[, 1:12]), gamma = 0.10)
png("output_plots_DTDnetwork/FemaleNetwork.png", width = 800, height = 800)
qgraph(glasso_full, layout = "spring", title = "Female Sample Network", theme = "colorblind", labels = colnames(group2_data[, 1:12]), vsize = 10, esize = 10)
dev.off()

# 运行NCT
print("startNCT!!!")
nct_test <- NCT(data1 = group1_data[, 1:12], 
                data2 = group2_data[, 1:12],
                gamma = 0.50, 
                test.edges = TRUE, 
                edges = "all",
                abs=TRUE,
                test.centrality = TRUE, 
                p.adjust.methods = "holm",
                paired = FALSE, 
                it = 1000)

summary(nct_test)

# 4. 计算中心性指标
# 使用 bootnet 包计算中心性指标
glasso_full <- abs(glasso_full)
centrality <- centralityPlot(glasso_full, include = c("Strength", "Closeness", "Betweenness"))

# 提取中心性指标
centrality_indices <- centrality_auto(glasso_full)
print(centrality_indices)
print('ok')

png("output_plots_DTDnetwork/centrality_indices.png", width = 800, height = 600)
print(centralityPlot(glasso_full, include = c("Strength", "Closeness", "Betweenness")))
dev.off()
print('good')

# 5. 桥接中心性分析
# 使用 networktools 包计算桥接中心性
bridge_results <- bridge(glasso_full)
summary(bridge_results)
print('yes')

# 可视化桥接中心性指标
png("output_plots_DTDnetwork/Bridge_Centrality1.png", width = 800, height = 600)
print(plot(bridge_results, include = c("Bridge Strength", "Bridge Betweenness", "Bridge Closeness")))
dev.off()
print('good')

# 6. 网络稳定性评估
# 使用 bootnet 包评估网络模型稳定性
boot_results <- bootnet(data, 
                        type = "person",
                        nBoots = 1000, # 引导抽样次数
                        default = "EBICglasso",
                        statistics = c("edge", "strength", "betweenness", "closeness", "distance")) # 使用 EBICglasso
print(boot_results)
# 检查稳定性结果
png("output_plots_DTDnetwork/BootResults.png", width = 800, height = 600)
print(plot(boot_results))
dev.off()
print('done')
