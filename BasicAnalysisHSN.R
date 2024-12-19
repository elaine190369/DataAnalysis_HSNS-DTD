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
data <- read.table("data.csv", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
data <- na.omit(data)
data <- data[, 1:10]
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

# 获取变量名
variable_names <- setdiff(names(data), "Total")

# 初始化存储结果的列表
t_test_results <- list()

# 对每一列（变量）进行 t 检验
for (var in variable_names) {
  # 提取变量在 high_group 和 low_group 中的值
  high_values <- high_group[[var]]
  low_values <- low_group[[var]]
  
  # 进行 t 检验
  t_test <- t.test(high_values, low_values)
  
  # 保存结果
  t_test_results[[var]] <- list(
    variable = var,
    p_value = t_test$p.value,
    statistic = t_test$statistic,
    conf_int = t_test$conf.int,
    mean_high = mean(high_values, na.rm = TRUE),
    mean_low = mean(low_values, na.rm = TRUE)
  )
}

# 打印结果
for (result in t_test_results) {
  cat("Variable:", result$variable, "\n")
  cat("  t-statistic:", result$statistic, "\n")
  cat("  p-value:", result$p_value, "\n")
  cat("  Confidence Interval:", result$conf_int, "\n")
  cat("  Mean (High Group):", result$mean_high, "\n")
  cat("  Mean (Low Group):", result$mean_low, "\n")
  cat("\n")
}

# # 计算每题的高低分组平均分差
# distinction <- sapply(1:(ncol(data) - 1), function(i) {
#   high_mean <- mean(high_group[[i]])
#   low_mean <- mean(low_group[[i]])
#   high_mean - low_mean
# })

# # 输出极端分组法的区分度
# names(distinction) <- colnames(data)[1:(ncol(data) - 1)]
# print("区分度:")
# print(distinction)

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

# 4. 效度分析-探索性因素分析（EFA）
# 4.1 判断是否适合因素分析：Barlett球形分析 + KMO法
kmo_result <- KMO(data)
print("\nKMO 检验结果：")
print(kmo_result)
bartlett_result <- cortest.bartlett(cor(data), n = nrow(data))
print("Bartlett 球形检验结果：")
print(bartlett_result)

# 4.2 主成分分析
data_scaled <- scale(data)
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
summary(pca_result)  # 显示主成分的方差比例（Proportion of Variance）和累计方差比例

# 查看主成分的载荷矩阵（即变量对各主成分的贡献）
print("===== 主成分载荷矩阵 =====")
print(pca_result$rotation)

# 查看样本的主成分得分（样本在主成分空间的投影）
print("===== 样本的主成分得分 =====")
print(head(pca_result$x))

# 提取每个主成分的方差解释比例
explained_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cumulative_variance <- cumsum(explained_variance)
print("===== 方差解释比例 =====")
print(data.frame(
  PC = paste0("PC", 1:length(explained_variance)),  
  Explained_Variance = round(explained_variance, 4),
  Cumulative_Variance = round(cumulative_variance, 4)
))

# 可视化
output_dir <- "output_plots_HSNbasic"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 绘制并保存图像
tryCatch({
  # 碎石图
  png(file.path(output_dir, "scree_plot.png"), width = 800, height = 600, res = 120)
  print(fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50)))
  dev.off()
  
  # 变量贡献图
  png(file.path(output_dir, "variable_contribution_plot.png"), width = 800, height = 600, res = 120)
  print(fviz_pca_var(pca_result, col.var = "contrib",
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE))
  dev.off()
  
 # 样本分布图
png(file.path(output_dir, "sample_distribution_plot.png"), width = 800, height = 600, res = 120)
print(fviz_pca_ind(pca_result, 
                   col.ind = "cos2", 
                   gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                   label = "none", 
                   geom.ind = "point",
                   repel = TRUE))
dev.off()
}, error = function(e) {
  message("Error in plotting: ", e$message)
})

# 4.3 平行分析
png(file.path(output_dir, "parallel_scree_plot.png"), width = 800, height = 600, res = 120)
fa_parallel <- fa.parallel(data, fa = "fa", n.iter = 100)
dev.off()
print("平行分析建议的因子数：")
print(fa_parallel)

# 因子数 = 4，进行探索性因子分析
efa_result <- fa(data, nfactors = 4, rotate = "oblimin", fm = "ml")
print("平行分析结果：")
print(efa_result)

# 打印因子载荷矩阵
print("因子载荷矩阵：")
print(efa_result$loadings)

#运行CFA
# 初始 CFA 模型
print('now start CFACFACFA')
cfa_model <- '
  # ML1: 因子1
  ML1 =~ V4 + V5 + V10

  # ML2: 因子2
  ML2 =~ V2 + V7

  # ML3: 因子3
  ML3 =~ V1 + V8

  # ML4: 因子4
  ML4 =~ V9
'

optimize_cfa <- function(model, data) {
  iteration <- 1  # 记录迭代次数
  repeat {
    # 拟合当前模型
    fit <- cfa(model, data = data)
    cat("\n--- Iteration", iteration, "Model Fit ---\n")
    fit_measures <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", "chisq", "df", "pvalue"))
    print(round(fit_measures, 3))
    
    std_loadings <- standardizedSolution(fit)
    loadings <- std_loadings[std_loadings$op == "=~", c("lhs", "rhs", "est.std")]
    
    # 打印当前所有变量的载荷系数
    cat("\n--- Factor Loadings ---\n")
    print(loadings)
    
    # 检查是否有载荷低于 0.4 的变量
    low_loading_vars <- loadings$rhs[loadings$est.std < 0.4]
    low_loading_values <- loadings$est.std[loadings$est.std < 0.4]
    
    # 如果没有低载荷变量，终止循环
    if (length(low_loading_vars) == 0) {
      cat("\nAll loadings are >= 0.4. Model optimization complete.\n")
      break
    }
    
    # 选择当前载荷最低的变量
    min_index <- which.min(low_loading_values)
    var_to_remove <- low_loading_vars[min_index]
    loading_to_remove <- low_loading_values[min_index]
    
    # 打印被剔除的变量及其载荷系数
    cat("\nVariable with the lowest loading (< 0.4) in Iteration", iteration, ":\n")
    cat("  - Variable:", var_to_remove, ", Loading:", round(loading_to_remove, 3), "\n")
    
    # 从模型中剔除该变量
    model <- gsub(paste0("\\+?\\s*", var_to_remove), "", model)  # 移除对应变量
    model <- gsub("=~\\s*\\+", "=~", model)                     # 移除多余的 "+"" 符号
    
    # 打印更新后的模型
    cat("\nUpdated Model after Iteration", iteration, ":\n", model, "\n")
    
    # 增加迭代计数
    iteration <- iteration + 1
  }
  
  # 返回最终拟合结果
  return(fit)
}

# 调用优化函数
final_fit <- optimize_cfa(cfa_model, data)

# 绘制最终模型路径图
semPaths(final_fit, "std", layout = "tree", residuals = FALSE, whatLabels = "std")

# 输出最终 CFA 模型的完整结果
cat("\n--- Final CFA Model Summary ---\n")
summary(final_fit, fit.measures = TRUE, standardized = TRUE)

# 提取并打印最终模型的标准化载荷
cat("\n--- Final Standardized Factor Loadings ---\n")
final_loadings <- standardizedSolution(final_fit)
final_loadings <- final_loadings[final_loadings$op == "=~", c("lhs", "rhs", "est.std")]
print(final_loadings)