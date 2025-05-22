# 加载必要的库
library(caret)
library(nnet)
library(ggplot2)
library(pROC)
library(doParallel)
library(VIM)    # 用于缺失值处理
library(Boruta) # 辅助特征选择

###############################
### 1. 数据加载与预处理函数 ###
###############################
load_and_preprocess <- function(file_path, log2_transform = TRUE) {
  # 读取数据
  raw_data <- read.csv(file_path, header = TRUE, check.names = FALSE, row.names = 1)
  
  # 转置为样本x基因格式
  data_matrix <- as.data.frame(t(raw_data))
  
  # 检查并处理缺失值（使用KNN插补）
  #if(any(is.na(data_matrix))) {
    #cat("检测到缺失值，正在进行KNN插补...\n")
    #data_matrix <- VIM::kNN(data_matrix, k = 5, imp_var = FALSE)
  #}
  
  # 自动判断是否需要log2转换
  if(log2_transform) {
    qx <- quantile(as.matrix(data_matrix), probs = c(0, 0.25, 0.5, 0.75, 0.99, 1))
    if(qx[5] > 100 || (qx[6] - qx[1] > 50 && qx[2] > 0)) {
      data_matrix <- log2(data_matrix + 1)
      cat("已应用log2(x+1)转换\n")
    } else {
      cat("无需log2转换\n")
    }
  }
  
  return(data_matrix)
}

##################################
### 2. 数据集划分函数（增强版） ###
##################################
create_train_test <- function(data_matrix, group_labels, train_ratio = 0.7, seed = 74521) {
  set.seed(seed)
  
  # 检查输入一致性
  if(nrow(data_matrix) != length(group_labels)) {
    stop("样本数与分组标签数量不匹配!")
  }
  
  # 创建分层抽样索引
  train_idx <- createDataPartition(
    y = group_labels,
    p = train_ratio,
    list = FALSE,
    times = 1
  )
  
  # 构建返回对象
  return(list(
    train_data = data_matrix[train_idx, ],
    test_data = data_matrix[-train_idx, ],
    train_labels = group_labels[train_idx],
    test_labels = group_labels[-train_idx]
  ))
}

#####################################
### 3. 神经网络模型核心函数（改进） ###
#####################################
run_nnet_analysis <- function(train_data, test_data, train_labels, test_labels,
                              output_dir = "NNET_Results",
                              num_genes = 30,
                              cv_folds = 5,
                              cv_repeats = 20,
                              n_cores = 4) {
  
  # 创建结果目录
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 设置并行计算
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl)) # 确保退出时释放资源
  
  ############################
  ### 第一阶段：特征选择 ###
  ############################
  
  # 使用Boruta算法辅助特征选择
  cat("正在进行特征初筛...\n")
  boruta_result <- Boruta(x = train_data, y = train_labels, doTrace = 0)
  important_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
  
  # 使用nnet进行精细特征排序
  ctrl <- trainControl(
    method = "repeatedcv",
    number = cv_folds,
    repeats = cv_repeats,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    allowParallel = TRUE
  )
  
  nnet_model <- train(
    x = train_data[, important_features],
    y = train_labels,
    method = "nnet",
    trControl = ctrl,
    tuneGrid = expand.grid(size = c(5, 10), decay = c(0.01, 0.1)),
    metric = "ROC",
    trace = FALSE
  )
  
  # 提取最终重要基因
  var_imp <- varImp(nnet_model)$importance
  selected_genes <- rownames(var_imp)[order(-var_imp$Overall)][1:num_genes]
  
  ############################
  ### 第二阶段：模型训练 ###
  ############################
  
  # 使用筛选后的特征集
  final_train <- train_data[, selected_genes]
  final_test <- test_data[, selected_genes]
  
  # 优化模型参数
  final_model <- train(
    x = final_train,
    y = train_labels,
    method = "nnet",
    trControl = trainControl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      savePredictions = "final",
      classProbs = TRUE,
      summaryFunction = twoClassSummary,
      allowParallel = TRUE
    ),
    tuneGrid = expand.grid(size = c(5, 10), decay = c(0.01, 0.1)),
    metric = "ROC",
    trace = FALSE
  )
  
  ############################
  ### 第三阶段：结果评估 ###
  ############################
  
  # 训练集交叉验证结果
  cv_perf <- final_model$results[
    final_model$results$size == final_model$bestTune$size &
      final_model$results$decay == final_model$bestTune$decay,
  ]
  
  # 测试集预测
  test_pred <- predict(final_model, newdata = final_test, type = "prob")
  test_roc <- roc(response = test_labels, predictor = test_pred$Case)
  
  # 可视化结果
  generate_plots(
    model = final_model,
    var_imp = var_imp,
    test_roc = test_roc,
    output_dir = output_dir
  )
  
  # 保存关键结果
  save_results(
    selected_genes = selected_genes,
    final_model = final_model,
    test_roc = test_roc,
    output_dir = output_dir
  )
  
  return(list(
    selected_genes = selected_genes,
    final_model = final_model,
    test_auc = test_roc$auc
  ))
}

################################
### 辅助函数：可视化与保存 ###
################################
generate_plots <- function(model, var_imp, test_roc, output_dir) {
  # 特征重要性图
  imp_plot <- ggplot(var_imp, aes(x = reorder(rownames(var_imp), Overall), y = Overall)) +
    geom_col(fill = "steelblue") +
    labs(title = "Feature Importance", x = "Gene", y = "Importance") +
    coord_flip() +
    theme_minimal()
  
  # ROC曲线图
  roc_plot <- ggroc(test_roc, color = "red") +
    geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
    labs(title = paste0("Test ROC (AUC = ", round(test_roc$auc, 3), ")"))
  
  # 保存图片
  ggsave(file.path(output_dir, "feature_importance.png"), imp_plot)
  ggsave(file.path(output_dir, "roc_curve.png"), roc_plot)
}

save_results <- function(selected_genes, final_model, test_roc, output_dir) {
  # 保存基因列表
  write.csv(selected_genes, file.path(output_dir, "selected_genes.csv"), row.names = FALSE)
  
  # 保存模型
  saveRDS(final_model, file.path(output_dir, "final_model.rds"))
  
  # 保存性能指标
  sink(file.path(output_dir, "performance.txt"))
  cat(paste("Test AUC:", round(test_roc$auc, 3), "\n"))
  print(final_model$results)
  sink()
}

##############################
### 4. 主执行流程 ###
##############################
# 参数设置
input_file <- "ht.csv"
log2_transform <- TRUE
train_ratio <- 0.7
num_genes <- 30
cv_folds <- 5
cv_repeats <- 20
n_cores <- 4

# 执行流程
tryCatch({
  # 数据加载与预处理
  cat("Step 1/4: 加载并预处理数据...\n")
  data_matrix <- load_and_preprocess(input_file, log2_transform)
  
  # 创建分组标签（假设前10例为Case，后10例为Control）
  group_labels <- factor(c(rep("Case", 10), rep("Control", 10)))
  
  # 划分训练测试集
  cat("Step 2/4: 划分数据集...\n")
  split_data <- create_train_test(data_matrix, group_labels, train_ratio)
  
  # 运行分析流程
  cat("Step 3/4: 运行神经网络分析...\n")
  results <- run_nnet_analysis(
    train_data = split_data$train_data,
    test_data = split_data$test_data,
    train_labels = split_data$train_labels,
    test_labels = split_data$test_labels,
    output_dir = "NNET_Results",
    num_genes = num_genes,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    n_cores = n_cores
  )
  
  # 最终输出
  cat("Step 4/4: 分析完成！结果已保存到NNET_Results目录\n")
  cat(paste("测试集AUC:", round(results$test_auc, 3), "\n"))
  
}, error = function(e) {
  cat("分析过程中出现错误:\n")
  cat(e$message, "\n")
  if(exists("cl")) stopCluster(cl) # 确保释放并行资源
})