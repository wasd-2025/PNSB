library(geneplotter)
library(copynumber)
#绘制代谢途径的热图
fitness123_annotated1 <- fitness123_annotated %>%
  mutate(time = recode(time, `0` = "T0-YE", `1` = "T1-Formate", `2` = "T2-Formate"))
heatmap_fitness <- function(data, key = TRUE, max_value = 6) {
  # 颜色设置
  heat_cols <- diverging_hcl(n = 7, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3))
  # 绘制热图
  levelplot(norm_gene_fitness_mean ~ factor(gene_name) * fct_rev(factor(time)),
            data = df_fdh,
            #par.settings = custom.colorblind(),colorkey = key,
            col.regions = colorRampPalette(heat_cols)(25),
            at = seq(-6, 6, 0.5), 
            aspect = "iso",
            xlab = "",
            ylab = "",
            scales = list(cex = 0.7, x = list(rot = 90)),
            panel = function(x, y, z, ...) {
              panel.levelplot(x, y, z, ...)
              panel.abline(v = seq_along(unique(x)) + 0.5, 
                           h = seq_along(unique(y)) + 0.5, 
                           col = "white", lwd = 2)
            }
  )
}
#基因位置图
genome_plot <- function(df, xlim = NULL, title = "", rot_labels = 0) {
  df <-  replace_na(df, list(strains_per_gene = 0))
  theme <- theme(
    axis.line = element_line(colour = grey(0.3)),
    axis.text = element_text(colour = grey(0.3), size = 0.7 * rel(1)), 
    axis.title = element_text(colour = grey(0.3)),
    plot.background = element_rect(fill = grey(0.7), colour = grey(0.85)),
    panel.background = element_rect(fill = grey(0.7), colour = grey(0.85)))
  
  if (!is.null(xlim))
    xscale = list(limits = xlim)
  else 
    xscale = list()
  xyplot(end/1000 ~ start/1000, df,
         groups = strand, cex = 0.7, lwd = 1,
         par.settings = theme, strains = df$strains_per_gene,
         scales = list(y = list(draw = FALSE), x = xscale),
         ylim = c(-3,2), xlab = "", ylab = "",
         gene_strand = df[["strand"]],
         gene_name = df[["gene_name"]],
         panel = function(x, y, strains, ...) {
           panel.geneplot(x, y, arrows = TRUE, tip = 0.1, rot_labels = rot_labels, ...)
           panel.text((x+y)/2, rep(0, length(x)), labels = strains, cex = 0.7)
           #panel.text(mean(xlim), 1.5, labels = title, col = 1)
         }
  )
}
df_fdh1 <- df_fdh
# 数据处理
df_fdh <- fitness123_annotated1 %>%
  dplyr::filter(locusId %in% group1_locus) %>%
  mutate(gene_name = gene)

# 检查并手动指定 gene_name
cat("正在为所有基因指定新的名字，请根据提示输入：\n")
for (i in seq_len(nrow(df_fdh))) {
  current_name <- ifelse("Name" %in% names(df_fdh), df_fdh$Name[i], NA)
  current_locusid <- df_fdh$locusId[i]
  prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
  
  # 确保输入不为空
  repeat {
    new_name <- readline(prompt = prompt_msg)
    if (!is.na(new_name) && new_name != "") {
      df_fdh$gene_name[i] <- new_name
      break
    } else {
      cat("输入不能为空，请重新输入：\n")
    }
  }
}

row_to_change <- 6
new_name <- "pufB2"
rc[row_to_change, "gene_name"] <- new_name
# 查看修改后的结果
print(rc)

# 绘制热图
plot_fdh_fit <- heatmap_fitness(df_fdh)

# 绘制基因组图
df_moco <- df_fdh %>%
  group_by(locusId, gene, scaffold, strand, start, end, gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")

plot_moco_g1 <- genome_plot(df_moco, xlim = c(4948.7, 4957.2), title = "chr 1")

# 打印图形到控制台
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)

# 导出热图到 SVG 文件
svg("heatmap_fdh_fit.svg", height = 6)
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
dev.off()

# 导出基因组图到 SVG 文件
svg("genome_plot_moco_g1_fdh.svg", height = 6)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()


# 数据处理\
cbb <- fitness123_annotated1 %>%
  dplyr::filter(locusId %in% group2_locus) %>%
  mutate(gene_name = gene)

# 检查并手动指定 gene_name
cat("正在为所有基因指定新的名字，请根据提示输入：\n")
for (i in seq_len(nrow(cbb))) {
  current_name <- ifelse("Name" %in% names(cbb), cbb$Name[i], NA)
  current_locusid <- cbb$locusId[i]
  prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
  
  # 确保输入不为空
  repeat {
    new_name <- readline(prompt = prompt_msg)
    if (!is.na(new_name) && new_name != "") {
      cbb$gene_name[i] <- new_name
      break
    } else {
      cat("输入不能为空，请重新输入：\n")
    }
  }
}
# 定义预设的基因集合
preset_genes <- c("rbcS", "rlp2")  # 替换为你的预设基因名称

df_fdh <- cbb
# 调整数据框的行顺序，将预设的基因集合放在最前面
df_fdh <- df_fdh %>%
  mutate(order = case_when(
    gene_name %in% preset_genes ~ 1,
    TRUE ~ 2
  )) %>%
  arrange(order) %>%
  select(-order)

# 绘制热图

plot_fdh_fit <- heatmap_fitness(df_fdh)

# 绘制基因组图
df_moco <- df_fdh %>%
  group_by(locusId, gene, scaffold, strand, start, end, gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")

plot_moco_g1 <- genome_plot(df_moco, xlim = c(1247, 1250.1), title = "chr 1")
plot_moco_g2 <- genome_plot(df_moco, xlim = c(4462.1, 4469.3), title = "chr 1")
plot_moco_g3 <- genome_plot(df_moco, xlim = c(3757.7, 3762.6), title = "chr 1")

# 打印图形到控制台
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
print(plot_moco_g3, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
# 导出热图到 SVG 文件
svg("heatmap_fdh_fitcbb.svg", height = 6)
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
dev.off()

# 导出基因组图1到 SVG 文件
svg("genome_plot_moco_g1cbb.svg", height = 4)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()

# 导出基因组图3到 SVG 文件
svg("genome_plot_moco_g3cbb.svg", height = 4)
print(plot_moco_g3, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()

# 导出基因组图2到 SVG 文件
svg("genome_plot_moco_g2cbb.svg", height = 4)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()



# ETC
etc <- fitness123_annotated1 %>%
  dplyr::filter(locusId %in% group3_locus) %>%
  mutate(gene_name = gene)

# 检查并手动指定 gene_name
cat("正在为所有基因指定新的名字，请根据提示输入：\n")
for (i in seq_len(nrow(etc))) {
  current_name <- ifelse("Name" %in% names(etc), etc$Name[i], NA)
  current_locusid <- etc$locusId[i]
  prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
  
  # 确保输入不为空
  repeat {
    new_name <- readline(prompt = prompt_msg)
    if (!is.na(new_name) && new_name != "") {
      etc$gene_name[i] <- new_name
      break
    } else {
      cat("输入不能为空，请重新输入：\n")
    }
  }
}
df_fdh <- etc
# 绘制热图

plot_fdh_fit <- heatmap_fitness(df_fdh)

# 绘制基因组图
df_moco <- df_fdh %>%
  group_by(locusId, gene, scaffold, strand, start, end, gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")

plot_moco_g1 <- genome_plot(df_moco, xlim =  c(4982.4, 4987.2), title = "chr 1")
plot_moco_g2 <- genome_plot(df_moco, xlim = c(4832.4, 4838), title = "chr 1")

# 打印图形到控制台
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)

# 导出热图到 SVG 文件
svg("heatmap_fdh_fitetc.svg", height = 6)
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
dev.off()

# 导出基因组图1到 SVG 文件
svg("genome_plot_moco_g1etc.svg", height = 4)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()

# 导出基因组图2到 SVG 文件
svg("genome_plot_moco_g2etc.svg", height = 4)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()

#ETC
df_ETC <- fitness123_annotated1 %>% 
  dplyr::filter(locusId %in% group3_locus) %>%
  mutate(gene_name = gene)
na_indices <- which(is.na(df_ETC$gene_name))
if (length(na_indices) > 0) {
  cat("以下行的 gene_name 为 NA，你可以为它们指定新的名字：\n")
  for (i in na_indices) {
    current_name <- ifelse("Name" %in% names(df_ETC), df_ETC$Name[i], NA)
    current_locusid <- df_ETC$locusId[i]
    prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
    new_name <- readline(prompt = prompt_msg)
    df_ETC$gene_name[i] <- new_name
  }
}
# 检查重复的 gene_name
duplicated_indices <- which(duplicated(df_ETC$gene_name) | duplicated(df_ETC$gene_name, fromLast = TRUE))
if (length(duplicated_indices) > 0) {
  cat("以下行的 gene_name 存在重复，你可以为它们指定新的名字：\n")
  for (i in duplicated_indices) {
    current_name <- df_ETC$gene_name[i]
    current_locusid <- df_ETC$locusId[i]
    # 找到重复的 gene_name 对应的所有 locusId
    duplicated_locusids <- df_ETC$locusId[df_ETC$gene_name == current_name]
    # 检查是否存在不同的 locusId
    if (length(unique(duplicated_locusids)) > 1) {
      # 提示用户输入新的 gene_name
      prompt_msg <- paste("第", i, "行的 gene_name（", current_name, "）存在重复，对应的 locusId 有：", paste(unique(duplicated_locusids), collapse = ", "), "，请为其指定新的名字（当前 locusId: ", current_locusid, "）：", sep = "")
      new_name <- readline(prompt = prompt_msg)
      df_ETC$gene_name[i] <- new_name
    }
  }
}
df_fdh <- df_ETC
plot_fdh_fit <- heatmap_fitness(df_fdh)
df_moco <- df_fdh %>% group_by(locusId, gene, scaffold, strand, start, end,gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")
# 分别绘制两组数据的图形
plot_moco_g1_group2 <- genome_plot(df_moco, xlim = c(4982.4, 4987.2), title = "chr 1")
plot_moco_g1_group1 <- genome_plot(df_moco, xlim = c(4832.4, 4838), title = "chr 1")
print(plot_fdh_fit, position = c(0,0.25,1,1.05), more = TRUE)
print(plot_moco_g1_group2,position = c(0.53,-0.1,0.8,0.45),more = TRUE)
print(plot_moco_g1_group1, position = c(0.2,-0.1,0.55,0.45), more = TRUE)






# RC
rc <- fitness123_annotated1 %>%
  dplyr::filter(locusId %in% group4_locus) %>%
  mutate(gene_name = gene)

# 检查并手动指定 gene_name
cat("正在为所有基因指定新的名字，请根据提示输入：\n")
for (i in seq_len(nrow(rc))) {
  current_name <- ifelse("Name" %in% names(rc), rc$Name[i], NA)
  current_locusid <- rc$locusId[i]
  prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
  
  # 确保输入不为空
  repeat {
    new_name <- readline(prompt = prompt_msg)
    if (!is.na(new_name) && new_name != "") {
      rc$gene_name[i] <- new_name
      break
    } else {
      cat("输入不能为空，请重新输入：\n")
    }
  }
}
df_fdh <- rc
# 绘制热图

plot_fdh_fit <- heatmap_fitness(df_fdh)

# 绘制基因组图
df_moco <- df_fdh %>%
  group_by(locusId, gene, scaffold, strand, start, end, gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")

plot_moco_g1 <- genome_plot(df_moco, xlim =  c(4006.1, 4009), title = "chr 1")
plot_moco_g2 <- genome_plot(df_moco, xlim = c(1230.1, 1237.8), title = "chr 1")

# 打印图形到控制台
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)

# 导出热图到 SVG 文件
svg("heatmap_fdh_fit_rc(1).svg", height = 6)
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
dev.off()

# 导出基因组图1到 SVG 文件
svg("genome_plot_moco_g1_rc(1).svg", height = 4)
print(plot_moco_g1, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()

# 导出基因组图2到 SVG 文件
svg("genome_plot_moco_g2_rc(1).svg", height = 4)
print(plot_moco_g2, position = c(0.2, -0.1, 0.8, 0.45), more = TRUE)
dev.off()
#RC（改版）
df_RC <- fitness123_annotated1 %>% 
  dplyr::filter(locusId %in% group4_locus) %>%
  mutate(gene_name = gene)
na_indices <- which(is.na(df_RC$gene_name))
if (length(na_indices) > 0) {
  cat("以下行的 gene_name 为 NA，你可以为它们指定新的名字：\n")
  for (i in na_indices) {
    current_name <- ifelse("Name" %in% names(df_RC), df_RC$Name[i], NA)
    current_locusid <- df_RC$locusId[i]
    prompt_msg <- paste("请为第", i, "行（Name: ", current_name, ", locusId: ", current_locusid, "）指定新的 gene_name：", sep = "")
    new_name <- readline(prompt = prompt_msg)
    df_RC$gene_name[i] <- new_name
  }
}
# 检查重复的 gene_name
duplicated_indices <- which(duplicated(df_RC$gene_name) | duplicated(df_RC$gene_name, fromLast = TRUE))
if (length(duplicated_indices) > 0) {
  cat("以下行的 gene_name 存在重复，你可以为它们指定新的名字：\n")
  for (i in duplicated_indices) {
    current_name <- df_RC$gene_name[i]
    current_locusid <- df_RC$locusId[i]
    # 找到重复的 gene_name 对应的所有 locusId
    duplicated_locusids <- df_RC$locusId[df_RC$gene_name == current_name]
    # 检查是否存在不同的 locusId
    if (length(unique(duplicated_locusids)) > 1) {
      # 提示用户输入新的 gene_name
      prompt_msg <- paste("第", i, "行的 gene_name（", current_name, "）存在重复，对应的 locusId 有：", paste(unique(duplicated_locusids), collapse = ", "), "，请为其指定新的名字（当前 locusId: ", current_locusid, "）：", sep = "")
      new_name <- readline(prompt = prompt_msg)
      df_RC$gene_name[i] <- new_name
    }
  }
}
df_fdh <- df_RC
plot_fdh_fit <- heatmap_fitness(df_fdh)
df_moco <- df_fdh %>% group_by(locusId, gene, scaffold, strand, start, end,gene_name) %>%
  summarize(strains_per_gene = min(unique(strains_per_gene)), .groups = "drop")
# 分别绘制两组数据的图形
plot_moco_g1_group2 <- genome_plot(df_moco, xlim = c(4006.1, 4009), title = "chr 1")
plot_moco_g1_group1 <- genome_plot(df_moco, xlim = c(1230.1, 1237.8), title = "chr 1")
print(plot_fdh_fit, position = c(0,0.25,1,1.05), more = TRUE)
print(plot_moco_g1_group2,position = c(0.53,-0.1,0.8,0.45),more = TRUE)
print(plot_moco_g1_group1, position = c(0.2,-0.1,0.55,0.45), more = TRUE)




#genename写回去
# 步骤1: 去重处理（保留每个locusid的第一条记录）
df_cbb_unique <- df_cbb[!duplicated(df_cbb$locusId), c("locusId", "gene_name")]
df_cbb_unique <- df_RC[!duplicated(df_RC$locusId), c("locusId", "gene_name")]
df_cbb_unique <- df_ETC[!duplicated(df_ETC$locusId), c("locusId", "gene_name")]
df_cbb_unique <- df_fdh[!duplicated(df_fdh$locusId), c("locusId", "gene_name")]
# 步骤2: 精确的列替换（不使用完整合并）
library(dplyr)

fitness123_annotated1 <- fitness123_annotated1 %>%
  # 只从df_cbb_unique获取gene_name列
  left_join(select(df_cbb_unique, locusId, gene_name), by = "locusId") %>%
  # 替换gene列
  mutate(gene = coalesce(gene_name, gene)) %>%
  # 移除临时列
  select(-gene_name)



df123 <- read.csv("123.csv", stringsAsFactors = FALSE)
df130 <- read.csv("130.csv", stringsAsFactors = FALSE)
#绘制代谢途径的热图
fitness123_annotated1 <- fitness123_annotated %>%
  mutate(time = recode(time, `0` = "0-6", `1` = "0-48", `2` = "6-48"))
heatmap_fitness <- function(data, key = TRUE, max_value = 8) {
  # 颜色设置
  heat_cols <- diverging_hcl(n = 8, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3))
  # 绘制热图
  levelplot(Log2FC_mean ~ factor(locusId) * fct_rev(factor(time)),
            data = df_fdh,
            #par.settings = custom.colorblind(),colorkey = key,
            col.regions = colorRampPalette(heat_cols)(25),
            at = seq(-5, 8, 0.6), 
            aspect = "iso",
            xlab = "",
            ylab = "",
            scales = list(cex = 0.7, x = list(rot = 90)),
            panel = function(x, y, z, ...) {
              panel.levelplot(x, y, z, ...)
              panel.abline(v = seq_along(unique(x)) + 0.5, 
                           h = seq_along(unique(y)) + 0.5, 
                           col = "white", lwd = 2)
            }
  )
}

df_fdh <- df130
# 绘制热图
plot_fdh_fit <- heatmap_fitness(df_fdh)

# 打印图形到控制台
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)


# 导出热图到 SVG 文件
svg("heatmap_fdh_fit130.svg", height = 6)
print(plot_fdh_fit, position = c(0, 0.25, 1, 1.05), more = TRUE)
dev.off()

