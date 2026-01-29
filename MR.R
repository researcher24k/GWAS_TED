setwd(dir="D:/mendelian randomization/test")  # 设置工作路径

# 加载需要的包
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(ieugwasr)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
# 读取单个暴露文件（改s为你的暴露文件名）
exp_dat <- fread("./BLV.txt")  # 确保列名正确
exp_dat <- as.data.frame(exp_dat)
head(exp_dat)
exp_data <- subset(exp_dat, p < 1e-5)
all_results=data.frame()
# 设置多个结局ID（示例）
outcome_ids <- paste0("ubm-b-", 1:3999)
outcome_ids <- c('ieu-b-2','ieu-b-7','ebi-a-GCST90038646','ieu-a-1025','ebi-a-GCST90038645','')
# 创建结果目录
dir.create("mendelian_AD", showWarnings = FALSE)
dir.create("temp", showWarnings = FALSE)
exp_data <- as.data.frame(exp_data)
# 处理暴露数据
exp_data <- format_data(
  exp_data,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  pval_col = "p",
  eaf_col = "freq",
  effect_allele_col = 'A1',
  other_allele_col = "A2"
)
library(rstatix)
trace(ld_clump_local, edit = T) # 修改后保存
setwd('../ukb/')
# 本地聚类
exp_clumped <- ld_clump_local(
  dplyr::tibble(rsid=exp_data$SNP, pval=exp_data$pval.exposure)[-1,],
  clump_kb = 1000,
  clump_r2 = 0.001,
  bfile='../ukb/maf/data_maf0.01_rs_ref',
  plink_bin='../ukb/plink/plink',
  clump_p = 1
)
library(TwoSampleMR)

exp_clumped <- ld_clump(
  dplyr::tibble(
    rsid = exp_data$SNP,
    pval = exp_data$pval.exposure
  ),
  clump_kb = 1000,
  clump_r2 = 0.001,
  pop = "EUR"
)

write.csv(exp_clumped,'exp_clumped.csv')
exp_clumped <- fread('./mendelian_BLINDANDVISIMPAIRMENT/exp_clumped_e-5.csv')
setwd('../GWAS_eye/mendelian_GRAVES_OPHT')
exp_data <- exp_data[exp_data$SNP %in% exp_clumped$rsid, ]
write.csv(exp_data,'exp_clumped_e-5.csv')

# 初始化结果数据框
all_results <- data.frame()
outcome_id <- outcome_ids[2]
access_token='eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMjczNzMyNzI5QHFxLmNvbSIsImlhdCI6MTc2ODkwNzMwMSwiZXhwIjoxNzcwMTE2OTAxfQ.tXVASawS4i0ApvU-D_fupCQl98fDGSwyjFCRydoqqYGbVKUmA9_xWZYDUHELV5uZWnQ_VMBHg9V7mLvp3TSOdxC1uIBqpFV-prXw0-54oMDUvjdqC2A-tA4u3Lxso60obWTC6Y0M4eFIqH69dwGoOJMBI_EiQvxbAX-uq0QaSednFtfSwCqE0lr0dPfP25bMXgj-1ydni3pfs-z9DQiq0D75xgPH1h_laZkNmQgHzXxxCSAvWnxyeCoqNzp0Sz584CiRsWVAF5PxEmJF6PDHbsFTdsxkIjZXnVmaGd-vuUMHD0oX-bE4s6rYweMgisCqubglBI6ydohnxcECoNTnLw'
ieugwasr::api_status()
# 载入包
library(dplyr) #这个包提供管道符号支持
library(ieugwasr)
usethis::edit_r_environ()

# 查看是否识别到token
ieugwasr::get_opengwas_jwt()
setwd('../')

# 保存最终结果
write.csv(all_results, "T1_volumn.csv", row.names = FALSE)
all_results <- T1_volumn
IVW_result <- all_results[all_results$method=='Inverse variance weighted',]
IVW_result <- IVW_result[IVW_result$pleiotropy>0.05,]

IVW_re <- IVW_result[IVW_result$analysis_direction%in%'反向',]
IVW_re <- IVW_re[IVW_re$pval<0.05,]
IVW_re$fdr <- p.adjust(IVW_re$pval,method='fdr')
IVW_result <- IVW_result[IVW_result$analysis_direction%in%'正向',]
IVW_result <- IVW_result[-which(IVW_result$id.outcome %in% IVW_re$id.exposure),]
IVW_result <- IVW_result[-which(IVW_result$id.outcome =='ubm-b-93'),]

IVW_result <- IVW_result[IVW_result$nsnp==10,]

IVW_result$fdr <- p.adjust(IVW_result$pval,method='fdr')
IVW_result_sig <- IVW_result[IVW_result$pval<0.05,]
IVW_result_sig_fdr<- IVW_result[IVW_result$fdr<0.1,]
IVW_result_sig_pos <- IVW_result_sig[IVW_result_sig$or>1,]
IVW_result_sig_neg <-  IVW_result_sig[IVW_result_sig$or<1,]
write.csv(IVW_result, "T1_volumn_new.csv", row.names = FALSE)

write.csv(IVW_result_sig, "T1_volumn_sig.csv", row.names = FALSE)


library(TwoSampleMR)
library(dplyr)

# 原正向分析结果中的显著结局 (假设已筛选出这些id)
sig_outcomes <- IVW_result_sig$id.outcome
all_reverse_results=data.frame()
# 反向MR分析函数
# 在循环开始前初始化结果存储
exp_data <- exp_dat
outcome_id <- outcome_ids[165]

N_AD_case    <- 2569	
N_AD_control <- 475933
N_AD_eff <- 4 / (1 / N_AD_case + 1 / N_AD_control)
exp_data$samplesize.exposure=N_AD_case+N_AD_control

#165'
#1-164 T1
#165-647 nonT1
#648-1049 cortical area
#1020-1325 cortical thickness

#1326-1366 intensity
#1367-1436 contraset
#1437
#1438-1451 T2*
#1452-1526 FA
#1527-1601 MO

#1602-1901 diffusivity
#1902-1976 ICVF
#1977-2051 OD
#2052-2126 ISOVF
#2127-2142 tfMRI
#2143-2218 node amplification

all_results <- list()
library(data.table)
library(dplyr)
library(stringr)

# 1. 读入所有结果文件
files <- c(
  "./cortical thickness_ivw_moderate.csv",
  "./cortical area_ivw_moderate.csv",
  "./T1_ivw_moderate.csv"
)

res_list <- lapply(files, fread)

# 2. 合并
res_all <- bind_rows(res_list)

# 3. 提取 outcome 中的 ubm 编号
# 假设 outcome 列叫 outcome / id.outcome / outcome_id（三种情况都兼容）
ubm_ids <- res_all %>%
  mutate(
    outcome_id = coalesce(
      outcome,
      id.outcome.x,
      outcome_id
    )
  ) %>%
  filter(str_detect(outcome_id, "^ubm-")) %>%
  distinct(outcome_id) %>%
  pull(outcome_id)

# 看一下
length(ubm_ids)
head(ubm_ids)

#1-164 T1
#165-647 nonT1
#648-1049 cortical area
#1020-1325 cortical thickness

#1326-1366 intensity
#1367-1436 contraset
#1437
#1438-1451 T2*
#1452-1526 FA
#1527-1601 MO

#1602-1901 diffusivity
#1902-1976 ICVF
#1977-2051 OD
#2052-2126 ISOVF
#2127-2142 tfMRI
#2143-2218 node amplification
#2219-3919 connectivity
all_results <- list()
library(TwoSampleMR)
library(MRPRESSO)
library(dplyr)
library(purrr)
library(data.table)
library(stringr)

# ----------------------------
# MR QC 过滤函数
# ----------------------------
mr_qc_filter <- function(dat,
                         F_cutoff = 10,
                         min_snps = 3,
                         exposure_type = c("binary", "quantitative")) {
  
  if (is.null(dat) || nrow(dat) < min_snps) return(NULL)
  if (any(is.na(dat$beta.exposure)) || any(is.na(dat$se.exposure))) return(NULL)
  
  exposure_type <- match.arg(exposure_type)
  
  if (exposure_type == "binary") {
    # Logistic GWAS: Wald-type F statistic
    dat$F <- (dat$beta.exposure^2) / (dat$se.exposure^2)
    dat$R2 <- NA_real_
    
  } else {
    # Quantitative exposure only
    if (any(is.na(dat$eaf.exposure)) || any(is.na(dat$samplesize.exposure)))
      return(NULL)
    
    N <- dat$samplesize.exposure[1]
    
    dat$R2 <- (2 * dat$eaf.exposure * (1 - dat$eaf.exposure) *
                 dat$beta.exposure^2) /
      (2 * dat$eaf.exposure * (1 - dat$eaf.exposure) *
         (dat$beta.exposure^2 + dat$se.exposure^2 * N))
    
    dat$F <- (N - 2) * dat$R2 / (1 - dat$R2)
  }
  
  dat <- subset(dat, F > F_cutoff)
  if (nrow(dat) < min_snps) return(NULL)
  
  return(dat)
}


# ----------------------------
# MR 核心分析
# ----------------------------
run_mr_core <- function(dat, exposure_id, outcome_id, direction) {
  
  methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
  res <- mr(dat, method_list = methods)
  res_or <- generate_odds_ratios(res)
  
  hetero <- mr_heterogeneity(dat)
  pleio  <- mr_pleiotropy_test(dat)
  steiger <- directionality_test(dat)
  
  # 合并结果
  final <- res_or %>%
    left_join(hetero, by = "method") %>%
    mutate(
      pleiotropy_p = pleio$pval,
      steiger_correct = steiger$correct_causal_direction[1],
      n_snps = length(unique(dat$SNP)),
      mean_F = mean(dat$F),
      exposure = exposure_id,
      outcome = outcome_id,
      direction = direction
    )
  
  return(final)
}
type_list <- list( T1 = 1:164, nonT1 = 165:647, cortical_area = 648:1019, cortical_thickness = 1020:1325, intensity = 1326:1366, contraset = 1367:1436, T2star = 1438:1451, FA = 1452:1526, MO = 1527:1601, diffusivity = 1602:1901, ICVF = 1902:1976, OD = 1977:2051, ISOVF = 2052:2126, tfMRI = 2127:2142, node_amplification = 2143:2218, connectivity = 2219:3919 )
# ----------------------------
# IDP 类型列表
# ----------------------------
fdr_by_type <- function(df, type_list) {
  
  library(dplyr)
  library(stringr)
  
  # ===============================
  # 1. 仅保留 IVW 结果用于 FDR
  # ===============================
  ivw_df <- df %>%
    filter(method == "Inverse variance weighted")
  
  # ===============================
  # 2. 根据方向选择用于分类的 ID
  #    并提取 ubm-b-XXXX 中的数字
  # ===============================
  ivw_df <- ivw_df %>%
    rowwise() %>%
    mutate(
      type_id_raw = case_when(
        direction == "BVL_to_AD" ~ id.outcome.x,
        direction == "AD_to_BVL" ~ id.exposure.x,
        TRUE ~ NA_character_
      ),
      type_id_num = as.numeric(str_extract(type_id_raw, "\\d+")),
      type = names(Filter(function(x) type_id_num %in% x, type_list))[1]
    ) %>%
    ungroup()
  
  # ===============================
  # 3. 按 type + direction 做 FDR
  # ===============================
  ivw_df <- ivw_df %>%
    group_by(type, direction) %>%
    mutate(fdr_pval = p.adjust(pval, method = "fdr")) %>%
    ungroup()
  
  # ===============================
  # 4. 回填 FDR 到所有 MR 方法
  # ===============================
  df <- df %>%
    left_join(
      ivw_df %>%
        select(id.exposure.x, id.outcome.x, direction, fdr_pval),
      by = c("id.exposure.x", "id.outcome.x", "direction")
    )
  
  # ===============================
  # 5. MR strength 分级
  # ===============================
  df <- df %>%
    mutate(
      mr_strength = case_when(
        !is.na(fdr_pval) &
          fdr_pval < 0.05 &
          pleiotropy_p > 0.05 &
       
          n_snps >= 3 ~ "strong",
        
        !is.na(fdr_pval) &
          fdr_pval < 0.05 &
          (pleiotropy_p <= 0.05 | !steiger_correct) &
          n_snps >= 3 ~ "moderate",
        
        TRUE ~ "weak"
      )
    )
  
  return(df)
}

# ----------------------------
# 分层 FDR 矫正函数

fdr_global <- function(df) {
  df <- df %>%
    mutate(fdr_pval = p.adjust(pval, method = "fdr")) %>%
    mutate(mr_strength = case_when(
      fdr_pval < 0.05 & pleiotropy_p > 0.05 & steiger_correct & n_snps >= 3 ~ "strong",
      fdr_pval < 0.05 & (pleiotropy_p <= 0.05 | !steiger_correct) & n_snps >= 3 ~ "moderate",
      TRUE ~ "weak"
    ))
  return(df)
}
# ----------------------------
# 主循环
# ----------------------------
all_results <- list()
idp_id <- outcome_ids[1]
sample_size <- 		63926
exp_ids <- c('ebi-a-GCST90018852','finn-b-H7_AMD','ebi-a-GCST90018814','finn-b-H7_BLINDANDVISIMPAIRMENT')
exp_id <- outcome_ids[1]
outcome_ids <- c('ebi-a-GCST90027158','ieu-b-7','ebi-a-GCST90038646','ieu-a-1025','ebi-a-GCST90038645','')

for (idp_id in outcome_ids[c(1:164, 648:1601)]) {
  
  message("Processing IDP: ", idp_id)
  
  tryCatch({
    # 正向 AD -> IDP
    exp_iv <- extract_instruments(outcomes = exp_id, p1 = 1e-5, clump = TRUE, r2 = 0.01, kb = 5000)
    
    idp_out <- extract_outcome_data(snps = exp_iv$SNP, outcomes = idp_id)
    dat_fwd <- harmonise_data(exp_iv, idp_out)
    dat_fwd <- subset(dat_fwd, mr_keep)
    dat_fwd <- mr_qc_filter(dat_fwd, exposure_type = "binary")
    if (!is.null(dat_fwd)) {
      res_fwd <- run_mr_core(dat_fwd, "glaucoma", idp_id, "glaucoma_to_	AD")
      all_results[[length(all_results)+1]] <- res_fwd
    }
    
    # 反向 IDP -> AD
    idp_iv <- extract_instruments(outcomes = idp_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
    if (!is.null(idp_iv) && nrow(idp_iv) >= 3) {
      ad_out <- format_data(exp_dat,
                            snps = idp_iv$SNP,
                            type = "outcome",
                            snp_col = "SNP",
                            beta_col = "b",
                            se_col = "se",
                            pval_col = "p",
                            eaf_col = "freq",
                            effect_allele_col = "A1",
                            other_allele_col = "A2")
      dat_rev <- harmonise_data(idp_iv, ad_out)
      dat_rev <- subset(dat_rev, mr_keep)
    
      dat_rev <- mr_qc_filter(dat_rev, exposure_type = "binary")
      if (!is.null(dat_rev)) {
        res_rev <- run_mr_core(dat_rev, idp_id, 'glAUCOMA', "	AD_to_glaucoma")
        all_results[[length(all_results)+1]] <- res_rev
      }
    }
    
  }, error = function(e) {
    message("Error at ", idp_id, ": ", e$message)
  })
}

# ----------------------------
# 合并结果 & 分层 FDR 矫正 & 分级
# ----------------------------
final_results <- bind_rows(all_results)
final_results <- fdr_by_type(final_results, type_list)

# 查看最终结果
head(final_results)

write.csv(
  final_results_df,
  file = "MR_AD_full_results_forward_reverse.csv",
  row.names = FALSE
)

all_results <- list()
for (idp_id in outcome_ids[c(1:164, 648:1601)]) {
  
  message("Processing IDP: ", idp_id)
  
  tryCatch({
    # 正向 AD -> IDP
    idp_out <- extract_outcome_data(snps = exp_data$SNP, outcomes = idp_id)
    dat_fwd <- harmonise_data(exp_data, idp_out)
    dat_fwd <- subset(dat_fwd, mr_keep)
    dat_fwd <- mr_qc_filter(dat_fwd, exposure_type = "binary",
                            N_case = N_AD_case, N_control = N_AD_control)
    if (!is.null(dat_fwd)) {
      res_fwd <- run_mr_core(dat_fwd, "AD", idp_id, "AD_to_IDP")
      all_results[[length(all_results)+1]] <- res_fwd
    }
    
    # 反向 IDP -> AD
    idp_iv <- extract_instruments(outcomes = idp_id, p1 = 1e-6, clump = TRUE, r2 = 0.001, kb = 10000)
    if (!is.null(idp_iv) && nrow(idp_iv) >= 3) {
      ad_out <- format_data(exp_dat,
                            snps = idp_iv$SNP,
                            type = "outcome",
                            snp_col = "SNP",
                            beta_col = "b",
                            se_col = "se",
                            pval_col = "p",
                            eaf_col = "freq",
                            effect_allele_col = "A1",
                            other_allele_col = "A2")
      dat_rev <- harmonise_data(idp_iv, ad_out)
      dat_rev <- subset(dat_rev, mr_keep)
      dat_rev$samplesize.exposure <- N_AD_case + N_AD_control
      dat_rev <- mr_qc_filter(dat_rev, exposure_type = "quantitative")
      if (!is.null(dat_rev)) {
        res_rev <- run_mr_core(dat_rev, idp_id, "AD", "IDP_to_AD")
        all_results[[length(all_results)+1]] <- res_rev
      }
    }
    
  }, error = function(e) {
    message("Error at ", idp_id, ": ", e$message)
  })
}

# ----------------------------
# 合并结果 & 全局 FDR 矫正 & 分级
# ----------------------------
final_results <- bind_rows(all_results)
final_results <- fdr_by_type(final_results,type_list)
write.csv(final_results,'BVL_IDP_e-5')
# 查看结果
head(final_results)

# 定义分组规则 (名称 = 起始编号:结束编号)
group_rules <- list(
  "T1" = 1:164,
  "nonT1" = 165:647,
  "cortical area" = 648:1049,
  "cortical thickness" = 1050:1325,
  "intensity" = 1326:1366,
  "contraset" = 1367:1436,
  "T2" = 1438:1451)

  "FA" = 1452:1526,
  "MO" = 1527:1601,
  "diffusivity" = 1602:1901,
  "ICVF" = 1902:1976,
  "OD" = 1977:2051,
  "ISOVF" = 2052:2126,
  "tfMRI" = 2127:2142,
  "node amplification" = 2143:2218,
  "connectivity" = 2219:3919
)

# 创建子表格列表
sub_tables <- list()
all_results <- final_results_df
# 提取分组数据
for (group_name in names(group_rules)) {
  # 获取当前组的编号范围
  num_range <- group_rules[[group_name]]
  
  # 生成对应的 id.outcome 格式 (ubm-b-编号)
  outcome_ids <- paste0("ubm-b-", num_range)
  
  # 从 all_result 中筛选匹配的行
  sub_df <- all_results[all_results$id.outcome.x %in% outcome_ids, ]
  
  # 添加分组名称列
  sub_df$group <- group_name
  
  # 存储到结果列表
  sub_tables[[group_name]] <- sub_df
  
  # 打印进度信息
  cat(sprintf("创建子表格: %-20s | 包含 %4d 行 | ID范围: ubm-b-%d 到 ubm-b-%d\n",
              group_name, nrow(sub_df), min(num_range), max(num_range)))
}

# 单独处理特殊编号 1437
special_id <- "ubm-b-1437"
special_df <- all_results[all_results$id.outcome == special_id, ]
if (nrow(special_df) > 0) {
  special_df$group <- "special_1437"
  sub_tables[["special_1437"]] <- special_df
  cat(sprintf("创建特殊子表格: special_1437    | 包含 %4d 行 | ID: %s\n", 
              nrow(special_df), special_id))
}

# 访问子表格示例:
# sub_tables[["T1"]]       # T1组数据
# sub_tables[["connectivity"]]  # connectivity组数据

# 2. 保存所有子表格为CSV文件
output_dir <- "mendelian_BLINDANDVISIMPAIRMENT"  # 指定输出目录
dir.create(output_dir, showWarnings = FALSE)  # 创建目录（如果不存在）

for (group_name in names(sub_tables)) {
  # 生成文件名
  file_name <- paste0(group_name, ".csv")
  file_path <- file.path(output_dir, file_name)
  
  # 保存CSV
  write.csv(sub_tables[[group_name]], file_path, row.names = FALSE)
  
  # 打印保存信息
  cat(sprintf("已保存: %-25s | 行数: %4d | 文件: %s\n",
              group_name, nrow(sub_tables[[group_name]]), file_path))
}

for (group_name in names(sub_tables)) {
  # 生成文件名
  group_name <- sub_tables[[group_name]]
  
}

# 保存最终结果
write.csv(all_results, "./mendelian_AMD/contraset.csv", row.names = FALSE)
# 创建目录存储结果

# 创建目录存储结果
result_dir <- "ivw_analysis_results"
dir.create(result_dir, showWarnings = FALSE)

# 初始化结果汇总表
summary_df <- data.frame(
  group = character(),
  total_IVW = integer(),
  pleiotropy_pass = integer(),
  reverse_sig = integer(),
  forward_total = integer(),
  forward_p05 = integer(),
  forward_fdr05 = integer(),
  forward_pos = integer(),
  forward_neg = integer(),
  stringsAsFactors = FALSE
)
for (group_name in names(sub_tables)) {
  
  group_data <- sub_tables[[group_name]]
  
  # 1. 仅提取 IVW
  IVW_result <- group_data[group_data$method == "Inverse variance weighted", ]
  
  if (nrow(IVW_result) == 0) next
  
  # 2. 标记 pleiotropy / heterogeneity（不删除）
  IVW_result$flag_no_pleio <- IVW_result$pleiotropy_p > 0.05
  IVW_result$flag_low_het  <- IVW_result$Q_pval > 0.05
  
  # 3. 处理反向 MR（仅标记）
  IVW_re <- IVW_result[
    IVW_result$analysis == "reverse" &
      IVW_result$pval < 0.05, 
  ]
  
  reverse_ids <- unique(IVW_re$id.exposure)
  
  IVW_result$reverse_flag <- IVW_result$id.outcome.x %in% reverse_ids
  
  # 4. 仅保留正向结果用于主分析
  IVW_forward <- IVW_result[IVW_result$analysis == 'forward', ]
  
  if (nrow(IVW_forward) == 0) next
  
  # 5. FDR（仅对正向）
  IVW_forward$fdr <- p.adjust(IVW_forward$pval, method = "fdr")
  
  # 6. 构建 Evidence Level（核心）
  IVW_forward$evidence_level <- with(
    IVW_forward,
    ifelse(
      pval < 0.05 & fdr < 0.05 &
        mean_F > 10 & nsnp >= 3 &
        flag_no_pleio & flag_low_het &
        !reverse_flag,
      "Strong",
      ifelse(
        pval < 0.05 & mean_F > 10,
        "Moderate",
        "Weak"
      )
    )
  )
  
  # 7. 提取不同证据等级
  res_strong   <- IVW_forward[IVW_forward$evidence_level == "Strong", ]
  res_moderate <- IVW_forward[IVW_forward$evidence_level == "Moderate", ]
  
  # 8. 保存结果
  write.csv(
    IVW_forward,
    file.path(result_dir, paste0(group_name, "_ivw_all_forward.csv")),
    row.names = FALSE
  )
  if(nrow(res_strong)>0){
  write.csv(
    res_strong,
    file.path(result_dir, paste0(group_name, "_ivw_strong.csv")),
    row.names = FALSE
  )}
  if(nrow(res_moderate>0)){
  write.csv(
    res_moderate,
    file.path(result_dir, paste0(group_name, "_ivw_moderate.csv")),
    row.names = FALSE
  )}
  
  # 9. 汇总统计
  group_summary <- data.frame(
    group = group_name,
    total_IVW = nrow(IVW_result),
    forward_total = nrow(IVW_forward),
    strong = nrow(res_strong),
    moderate = nrow(res_moderate),
    weak = sum(IVW_forward$evidence_level == "Weak"),
    reverse_flag = sum(IVW_forward$reverse_flag)
  )
  
  summary_df <- rbind(summary_df, group_summary)
  
  cat(sprintf(
    "分组: %-20s | Strong: %d | Moderate: %d | Reverse_flag: %d\n",
    group_name,
    nrow(res_strong),
    nrow(res_moderate),
    sum(IVW_forward$reverse_flag)
  ))
  
  rm(IVW_result, IVW_forward, IVW_re)
}

# 保存与打印汇总
write.csv(summary_df, file.path(result_dir, "analysis_summary.csv"), row.names = FALSE)

cat("\n===== 分析结果汇总 =====\n")
print(summary_df)

library(forestplot)
library(dplyr)

# 模拟创建与您数据结构一致的示例数据
mr_results <- IVW_result_sig_fdr


# 自定义主题美化函数
compact_theme <- function() {
  fpTxtGp(
    label = gpar(fontfamily = "sans", cex = 0.85),
    ticks = gpar(cex = 0.8, col = "gray30"),
    xlab = gpar(cex = 0.9, fontface = "bold"),
    title = gpar(cex = 1.1, fontface = "bold", lineheight = 0.8)
  )
}

# 数据预处理管道
formatted_data <- mr_results %>%
  arrange(or) %>%
  mutate(
    or_ci = sprintf("%.2f (%.2f–%.2f)", or, or_lci95, or_uci95),
    p_formatted = case_when(
      pval < 0.001 ~ format(pval, scientific = TRUE, digits = 1),
      pval < 0.01 ~ sprintf("%.3f", pval),
      TRUE ~ sprintf("%.2f", pval)
    ),
    outcome = factor(outcome, levels = outcome)  # 锁定排序
  ) 
formatted_data <- mr_results %>%
  arrange(or) %>%
  mutate(
    or_ci = sprintf("%.2f (%.2f–%.2f)", or, or_lci95, or_uci95),
    p_formatted = case_when(
      pval < 0.001 ~ format(fdr, scientific = TRUE, digits = 1),
      pval < 0.01 ~ sprintf("%.3f", fdr),
      TRUE ~ sprintf("%.2f", fdr)
    ),
    outcome = factor(outcome, levels = outcome)  # 锁定排序
  ) 

# 构建紧凑标签矩阵
label_matrix <- cbind(
  c("Outcome", as.character(formatted_data$outcome)),
  c("OR (95% CI)", formatted_data$or_ci),
  c("P Value", formatted_data$p_formatted)
)
total_rows <- 1 + nrow(formatted_data)

# 生成森林图
forestplot(
  labeltext = label_matrix,
  mean = c(NA, formatted_data$or),
  lower = c(NA, formatted_data$or_lci95),
  upper = c(NA, formatted_data$or_uci95),
  is.summary = c(TRUE, rep(FALSE, nrow(formatted_data))),
  graph.pos = 2,  # 将图形移到中间列
  xlog = TRUE,
  col = fpColors(
    box = "#2C77BF",
    line = "#619CFF",
    summary = "firebrick",
    hrz_lines = "#4D4D4D"
  ),
  boxsize = 0.18,  # 更小的估计点
  lwd.ci = 1.8,    # 更细的置信区间线
  lineheight = unit(6, "mm"),  # 压缩行高
  txt_gp = compact_theme(),
  grid = structure(c(0.95, 1, 1.05), gp = gpar(lty = 3, col = "#F0F0F0")),
  xlab = "Odds Ratio(log)",
  title = expression(bold("MR Analysis: Exposure Effects on Outcomes")),
  
  hrzl_lines = list(
    "1" = gpar(lwd = 1.2, col = "#404040"),
    "2" = gpar(lwd = 0.6, col = "#A0A0A0"),
    '7' = gpar(lwd = 1.2, col = "#404040")  # 最后一行下方
  ),
  ci.vertices = TRUE,
  ci.vertices.height = 0.12,
  clip = c(0.9, 1.1),  # 紧凑的坐标范围
  mar = unit(c(2, 1, 1, 1), "mm")  # 最小化边距
)

# 添加副标题（需要配合grid包）
grid.text("Reference line: OR = 1", 
          x = unit(0.78, "npc"), y = unit(0.03, "npc"),
          gp = gpar(fontsize = 8, col = "gray40"))

library(meta)
library(metasens)


library(tidyverse)
library(forestplot)

r <- function(file_path) {
  # 文件处理管道
  combined_df <- list.files(path = './', pattern = "\\.csv$", full.names = TRUE) %>% 
    map_dfr(~{
      read_csv(.x) %>% 
        mutate(source = str_remove(basename(.x), "\\.csv$")) %>%  # 移除扩展名
        #filter(method == 'Inverse variance weighted') %>%
        #filter(pleiotropy > 0.05) %>%
        #filter(analysis_direction == '正向') %>%
        mutate(fdr = p.adjust(pval, method = 'fdr')) #%>%
      # filter(fdr < 0.05)
    })
  outcome_info <- combined_df %>%
    mutate(IDP_ID = as.integer(str_extract(id.outcome, "[0-9]+$")))
  # 把 IDP_short_name 加入到 outcome_info 中
  merged_info <- outcome_info %>%
    left_join(IDP, by =  "IDP_ID")
  write.csv(merged_info,'all.csv')
  # 生成带分组标题的数据结构
  formatted_data <- combined_df %>%
    group_by(source) %>%
    group_map(~{
      # 添加分组标题行
      header <- tibble(
        outcome = unique(.y$source),
        or = NA_real_, or_lci95 = NA_real_, or_uci95 = NA_real_,
        fdr = NA_real_, or_ci = "", p_formatted = ""
      )
      # 处理数据行
      body <- .x %>%
        arrange(or) %>%
        mutate(
          or_ci = sprintf("%.2f (%.2f–%.2f)", or, or_lci95, or_uci95),
          p_formatted = ifelse(fdr < 0.001, 
                               format(fdr, scientific = TRUE, digits =1),
                               sprintf("%.3f", fdr))
        )
      bind_rows(header, body)
    }) %>% 
    bind_rows() %>%
    mutate(outcome = factor(outcome, levels = unique(outcome)))
  
  # 构建可视化元素
  label_matrix <- cbind(
    c("Study", as.character(formatted_data$outcome)),
    c("OR (95% CI)", formatted_data$or_ci),
    c("FDR Value", formatted_data$p_formatted)
  )
  
  # 动态计算分割线位置
  split_lines <- c(1, which(is.na(formatted_data$or)), nrow(formatted_data)+1)
  
  # 增强型绘图主题
  enhanced_theme <- function() {
    fpTxtGp(
      label = gpar(fontfamily = "sans", col = c("black", rep("gray20", nrow(formatted_data)))),
      summary = gpar(fontface = "bold", col="navy"),
      xlab = gpar(cex = 0.9),
      title = gpar(cex = 1.2, fontface = "bold")
    )
  }
  
  # 生成森林图
  forestplot(
    labeltext = label_matrix,
    mean = c(NA, formatted_data$or),
    lower = c(NA, formatted_data$or_lci95),
    upper = c(NA, formatted_data$or_uci95),
    is.summary = c(TRUE, is.na(formatted_data$or)),
    graph.pos = 2,
    xlog = TRUE,
    col = fpColors(
      box = c(rep("#2196F3", nrow(formatted_data))),
      summary = "navy",
      hrz_lines = "#37474F"
    ),
    boxsize = 0.2,
    lwd.ci = 1.5,
    lineheight = unit(7, "mm"),
    txt_gp = enhanced_theme(),
    grid = c(1, 1.05),
    xlab = "Odds Ratio (log scale)",
    title = "Meta-analysis of MR Results by Data Source",
    #hrzl_lines = list(gpar(lwd=1, col="gray50")),
    
    ci.vertices = TRUE,
    ci.vertices.height = 0.15,
    clip = range(na.omit(c(formatted_data$or_lci95, formatted_data$or_uci95))),
    mar = unit(c(3, 1, 2, 1), "mm")
  )
  
  return(combined_df)
}

library(ggplot2)
library(dplyr)
library(ggrepel)

# 假设你的数据叫 result_df
# 如果你还没有添加 FDR，可以先计算
result_df <- merged_info %>%
  mutate(
    log_OR = log(or)
    #neg_log10_FDR = -log10(p.adjust(pval, method = "fdr"))
  )%>%filter(method == 'Inverse variance weighted') %>%
  filter(pleiotropy > 0.05)

# 选出最显著的前 5 个点用于标注
top_hits <- result_df%>%
  filter(fdr < 0.05)
#%>% arrange(fdr)%>%   slice_head(n = 5)

# 绘图
ggplot(result_df, aes(x = log_OR, y = -log10(fdr), color = Category_name)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(
    data = top_hits,
    aes(label = IDP_short_name),
    size = 3.5,
    max.overlaps = Inf
  ) +
  labs(
    x = "ln(OR)",
    y = "-log10(FDR)",
    title = "Volcano Plot of IDPs",
    color = "Category"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )



head(result_df)
merged_info <- all
data_wide <- merged_info%>%
  pivot_wider(names_from = method, values_from = c(pval, b))
data_wide = data_wide[,c("IDP_short_name","b_Inverse variance weighted",
                         "b_MR Egger","b_Weighted median","b_Weighted mode","b_Simple mode","pval_Inverse variance weighted",'Category_name','fdr')]
names(data_wide) = c("IDP_trait","IVW_beta","MR_Egger_beta","Weighted_median_beta","Weighted_mode_beta","Simple_mode_beta","IVW_pval",'Category','fdr')
library(dplyr)

df <- data_wide %>% 
  group_by(IDP_trait) %>% 
  summarise(
    IVW_beta = first(na.omit(IVW_beta)),
    MR_Egger_beta = first(na.omit(MR_Egger_beta)),
    Weighted_median_beta = first(na.omit(Weighted_median_beta)),
    Weighted_mode_beta = first(na.omit(Weighted_mode_beta)),
    Simple_mode_beta = first(na.omit(Simple_mode_beta)),
    IVW_pval = first(na.omit(IVW_pval)),
    Category=first(na.omit(Category)),
    fdr=min(fdr)
  )
data = df[df$fdr<0.05,]

data = data[order(data$IVW_beta,decreasing = T),]
hhh=data %>% as.data.frame()#数据转化为数据框
##设置工作路径
# setwd("/Users/ww/Desktop/R code/Figure/环状图")
#将数据的第一列转变为行名
rownames(hhh)=hhh[,1]
#删除第一列数据
hhh=hhh[,-1]
hhh=na.omit(hhh)
library(circlize)
#对不同表型设置颜色盘，前面6个是相关结果对应的颜色范围，最后一个是对class的颜色定义。
#需要绘制n个结果就需要输入n+1个子list，多出来的一个是class的颜色定义
col_fun = list(col_1 = colorRamp2(c(-0.1, 0, 0.1), 
                                  c("#0775b4", "white", "#a54b05")),
               col_2 = colorRamp2(c(-0.2, 0, 0.2),
                                  c( "#41ab5d", "white","#f25757")),
               col_3 = colorRamp2(c(-0.1,0,0.1),
                                  c( "#4F9AE9", "white", "#E9499F")),
               col_4 = colorRamp2(c(-0.1,0,0.1),
                                  c( "#181818","#181815", "white")),
               col_5 = colorRamp2(c(-0.1,0,0.1),
                                  c("#dbaf70", "white", "#dbaf75")),
               col_6 = colorRamp2(c(0.01, 0.05, 1),
                                  c( "#08306b","white","#08304b"))
)

#绘图

pdf("grave_circle_new.pdf", height = 12, width = 12) #设置输出的pdf文件长宽大小
circos.par$gap.degree <- 70 #相邻扇区的角度度数
circos.par$start.degree <- 20 #扇区从圆的顶部开始的几度角度处开始
circos.par$track.margin <- c(0.001, 0.001) #设置环形图中轨道之间的边距（间距）
for (i in 1:6) { #有几个结果就需要把i循环的上限改成相应的数字加1
  data_tmp <- as.matrix(hhh[,i]) #提取相应的结果
  if (i == 1) {
    rownames(data_tmp) <- rownames(hhh) #第一圈外围的图例
  }
  colnames(data_tmp) <- colnames(hhh)[i] #列名
  
  if (i <= 6) {
    #开始绘制圆形热图
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]], #提取上述设定的颜色
                   rownames.side = "outside", #设定名字在扇形的外面
                   cluster = T, #分类和分组，能够让后面的class与相应的traits对应
                   cell.border = "white", #单元格之间的边界颜色
                   track.height = 0.06) #热图的轨道高度
  } else {
    circos.heatmap(data_tmp, 
                   col = col_fun[[i]],
                   rownames.side = "outside", 
                   cluster = F,
                   track.height = 0.04) 
  }
}
library(ComplexHeatmap)#设置图例
lgd1 <- Legend(title = colnames(hhh)[1], #相应结果的标题名字
               border = "black", #设置图例边框的颜色为黑色
               grid_height = unit(3, "mm"), #每个图例项之间的垂直间距为 3 毫米
               legend_width = unit(20, "mm"), #参数设置图例的宽度为 20 毫米
               at = c(-0.5, 0.5),  #显示图例中各个cutoff的数值
               title_position = "topcenter", #标题的位置
               col_fun = col_fun[[1]], #颜色的图例
               direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold")) #图例的方向

lgd2 <- Legend(title = colnames(hhh)[2], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at = c(-0.5, 0.5), 
               title_position = "topcenter",
               col_fun = col_fun[[2]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd3 <- Legend(title = colnames(hhh)[3], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at =c(-0.05, 0.05),  
               title_position = "topcenter",
               col_fun = col_fun[[3]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd4 <- Legend(title = colnames(hhh)[4], border = "black", 
               grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at = c(-0.05, 0.05),
               title_position = "topcenter",
               col_fun = col_fun[[4]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd5 <- Legend(title = colnames(hhh)[5], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at =c(0.05,0.1), 
               title_position = "topcenter",
               col_fun = col_fun[[5]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

lgd6 <- Legend(title = colnames(hhh)[6], border = "black", grid_height = unit(3, "mm"),
               legend_width = unit(20, "mm"),
               at = c(0.05,0.1), 
               title_position = "topcenter",
               col_fun = col_fun[[6]], direction = "horizontal",title_gp = gpar(fontsize = 8, fontface = "bold"))

pd <- packLegend(lgd1, lgd2, lgd3, lgd4, lgd5, lgd6, row_gap = unit(1, "mm")) #对上述图例进行整合成一个图例
draw(pd, x = unit(0.6, "npc"), y = unit(0.75, "npc")) #设置图例所在位置

dev.off() #导出pdf
circos.clear() #清除工作环境中已有的环图信息

plot_df <- data.frame(
  region = c("insula", "precentral", "postcentral", "thalamus"),
  layer  = c("T1", "non-T1", "cortical", "subcortical"),
  evidence = c("strong", "moderate", "strong", "moderate")
)

library(ggseg3d)
library(ggseg)
data(dk_3d)
ggseg3d(
  atlas = dk_3d,        # Desikan atlas
  data = plot_df,
  region = "region",
 # colour = "layer",
  opacity = "evidence"
)

plot_df$alpha <- ifelse(plot_df$evidence == "strong", 1, 0.6)

ggseg3d(
  atlas = dk_3d,
  data = plot_df,
  region = "region",
  colour = "layer",
  opacity = "alpha"
)


