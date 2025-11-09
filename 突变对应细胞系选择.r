#!/usr/bin/env Rscript
# ============================================================
# 最小细胞系覆盖算法（ILP精确解版）
# 输入: mutations_subsetted_NAsdropped.csv
# 输出: optimized_cover_cells_ILP.csv 等三份结果
# ============================================================

# ---------- 0. 环境准备 ----------
needed <- c("dplyr", "tibble", "purrr", "ompr", "ompr.roi", "ROI.plugin.glpk")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(dplyr)
library(tibble)
library(purrr)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)

message("✅ 所需R包已加载。")

# ---------- 1. 输入文件 ----------
file_path <- "mutations_subsetted_NAsdropped.csv"
if (!file.exists(file_path)) stop("❌ 文件不存在：", file_path)
df <- tryCatch(
  read.csv(file_path, header = TRUE, stringsAsFactors = FALSE),
  error = function(e) read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
)
message("✅ 已读入数据：", file_path, "（行=", nrow(df), ", 列=", ncol(df), "）")

# ---------- 2. 识别基因列 ----------
preferred_cols <- c("gene", "Hugo_Symbol", "hgnc_name", "hgnc", "gene_symbol", "symbol")
found_cols <- intersect(preferred_cols, names(df))
if (length(found_cols) == 0) stop("❌ 未找到 gene/Hugo_Symbol/hgnc 等列。")

gene_col <- if ("gene" %in% names(df)) "gene" else found_cols[1]
df <- df %>% mutate(.gene_sym = toupper(as.character(.data[[gene_col]])))
message("✅ 使用列 '", gene_col, "' 作为基因列。")

# ---------- 3. 目标基因 ----------
target_genes <- c(
  "AKT1","ATM","BARD1","BRAF","BRCA1","BRCA2","BRIP1","CDH1","CDK4","CDKN2A",
  "CHEK2","EGFR","ERBB2","ESR1","FANCC","MRE11","NBN","NF1","NOTCH1","NTHL1",
  "NTRK1","NTRK2","NTRK3","PALB2","PIK3CA","PTEN","RAD50","RAD51B","RAD51C",
  "RAD51D","RECQL4","RET","RINT1","SLX4","SMARCA4","STK11","TP53","XRCC2","ZNF331"
)
target_genes_u <- unique(toupper(target_genes))

# ---------- 4. 过滤并构建 cell_gene_list ----------
if (!"cell_line_display_name" %in% names(df)) {
  stop("❌ 未找到细胞系列 (cell_line_display_name)。请检查输入文件。")
}

df_target <- df %>% filter(.gene_sym %in% target_genes_u)
message("✅ 匹配到目标基因的行数：", nrow(df_target))

cell_gene_list <- df_target %>%
  group_by(cell_line_display_name) %>%
  summarise(genes = list(unique(.gene_sym)), .groups = "drop")

cell_gene_list <- cell_gene_list %>%
  filter(length(intersect(genes[[1]], target_genes_u)) > 0)
message("✅ 构建了 cell_gene_list，包含细胞系数：", nrow(cell_gene_list))

# ---------- 5. 构建覆盖矩阵 ----------
genes <- unique(target_genes_u)
cells <- cell_gene_list$cell_line_display_name
A <- sapply(genes, function(g)
  sapply(cell_gene_list$genes, function(gs) g %in% gs)
)
A <- t(A) # 行=基因, 列=细胞系
rownames(A) <- genes
colnames(A) <- cells

message("✅ 构建覆盖矩阵完成。维度：", nrow(A), " 基因 × ", ncol(A), " 细胞系。")

# ===============================================
# 🔧 从内存矩阵 A 求解 ILP 模型（最终版）
# ===============================================

suppressPackageStartupMessages({
  library(ompr)
  library(ompr.roi)
  library(ROI.plugin.glpk)
  library(dplyr)
})

message("⚙️ 使用内存矩阵 A 求解最小细胞系集合...")

n_cells <- ncol(A)
n_genes <- nrow(A)
cells <- colnames(A)
genes <- rownames(A)

model <- MIPModel() %>%
  add_variable(x[i], i = 1:n_cells, type = "binary") %>%
  set_objective(sum_expr(x[i], i = 1:n_cells), "min")

# 每个基因至少由一个细胞系覆盖
for (j in seq_len(n_genes)) {
  coeffs <- as.numeric(A[j, ])
  if (all(coeffs == 0)) {
    stop(sprintf("❌ 基因 %s 无法被任何细胞系覆盖。", genes[j]))
  }
  idx_nonzero <- which(coeffs != 0)
  model <- model %>%
    add_constraint(sum_expr(coeffs[idx_nonzero[k]] * x[idx_nonzero[k]], k = seq_along(idx_nonzero)) >= 1)
}

message("✅ 模型构建完成，开始求解...")

result <- solve_model(model, with_ROI(solver = "glpk", verbose = TRUE))
status <- result$status
message(sprintf("📊 求解状态: %s", status))

# 输出结果
solution <- get_solution(result, x[i]) %>%
  filter(value > 0.5) %>%
  mutate(cell_line = cells[i]) %>%
  select(cell_line)

out_file <- "ILP_selected_cells.csv"
write.csv(solution, out_file, row.names = FALSE)

message(sprintf("✅ 已保存结果至: %s", out_file))
message(sprintf("📉 最少细胞系数目: %d", nrow(solution)))
