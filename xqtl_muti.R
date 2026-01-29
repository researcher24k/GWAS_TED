#!/usr/bin/env Rscript

# ============================================================================
# 共定位分析脚本 - 自动检测基因组版本 - 跳过已有结果 - 并行版本
# ============================================================================

# 解决包版本冲突
.libPaths(c(
  Sys.getenv("R_LIBS_USER"),
  .libPaths()
))

# 确保使用正确版本的Bioconductor包
if (packageVersion("BiocGenerics") < "0.51.2") {
  message("WARNING: BiocGenerics version too old, attempting to update...")
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("BiocGenerics", update = TRUE, ask = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(xQTLbiolinks)
  library(dplyr)
  library(optparse)
  library(xQTLbiolinks)
  library(coloc)
  library(hyprcoloc)
  library(data.table)
  library(stringr)
  library(R.utils)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(VariantAnnotation)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# ============================================================================
# 命令行参数解析
# ============================================================================
option_list <- list(
  make_option(c("-i", "--input_dir"), type="character", default="./vcf/converted_gwas",
              help="Directory containing GWAS files (*.smr.txt.gz) [default: %default]"),
  make_option(c("-o", "--output_dir"), type="character", default="./coloc_result",
              help="Output directory [default: %default]"),
  make_option(c("-p", "--pvalue"), type="numeric", default=5e-6,
              help="P-value threshold for sentinel SNPs [default: %default]"),
  make_option(c("-f", "--force"), action="store_true", default=FALSE,
              help="Force rerun even if results exist [default: %default]"),
  make_option(c("-n", "--ncores"), type="integer", default=6,
              help="Number of parallel cores [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 创建输出目录
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# ============================================================================
# 基因组版本检测函数(从Python移植)
# ============================================================================
detect_genome_version <- function(gwas_df) {
  # hg19/GRCh37 锚点 SNP
  anchors <- list(
    'rs3131972' = 752721,
    'rs530212009' = 628154,
    'rs13303240' = 928720,
    'rs12124819' = 1022045,
    'rs41285790' = 2507108,
    'rs12082461' = 10007631,
    'rs12119391' = 20002875,
    'rs2236357' = 50002830,
    'rs12039209' = 100000350,
    'rs10915170' = 200000302
  )
  
  # 检查是否有位置信息
  if (!all(c("SNP", "pos", "chr") %in% colnames(gwas_df))) {
    warning("Missing position columns, assuming GRCh38")
    return(list(version = "grch37", grch37To38 = TRUE))
  }
  
  # 取前100000行进行检测
  sample_size <- min(nrow(gwas_df), 100000)
  df_sample <- head(gwas_df, sample_size)
  
  # 找到匹配的锚点SNP
  matched <- df_sample[df_sample$SNP %in% names(anchors), ]
  
  if (nrow(matched) == 0) {
    message("No anchor SNPs found, assuming GRCh38")
    return(list(version = "grch37", grch37To38 = TRUE))
  }
  
  # 检查位置匹配度
  correct_count <- 0
  for (i in 1:nrow(matched)) {
    snp_id <- matched$SNP[i]
    actual_pos <- as.numeric(matched$pos[i])
    expected_pos <- anchors[[snp_id]]
    
    if (!is.na(actual_pos) && abs(actual_pos - expected_pos) < 5) {
      correct_count <- correct_count + 1
    }
  }
  
  # 如果有2个以上锚点匹配,判定为GRCh37
  is_grch37 <- correct_count >= 2
  
  if (is_grch37) {
    message("Detected genome version: GRCh37/hg19 (", correct_count, " anchors matched)")
    message("Will convert to GRCh38 for analysis")
    return(list(version = "grch37", grch37To38 = TRUE))
  } else {
    message("Detected genome version: GRCh38/hg38")
    return(list(version = "grch38", grch37To38 = FALSE))
  }
}

# ============================================================================
# P值清理函数
# ============================================================================
clean_pval <- function(df, p_col) {
  df <- df[!is.na(df[[p_col]]), ]
  df[[p_col]][df[[p_col]] <= 0] <- 1e-300
  df[[p_col]][df[[p_col]] > 1] <- 1
  return(df)
}

# ============================================================================
# GTEx 组织列表
# ============================================================================
tissues <- c(
  "Brain - Amygdala",
  "Brain - Anterior cingulate cortex (BA24)",
  "Brain - Caudate (basal ganglia)",
  "Brain - Cerebellar Hemisphere",
  "Brain - Cerebellum",
  "Brain - Cortex",
  "Brain - Frontal Cortex (BA9)",
  "Brain - Hippocampus",
  "Brain - Hypothalamus",
  "Brain - Nucleus accumbens (basal ganglia)",
  "Brain - Putamen (basal ganglia)",
  "Brain - Spinal cord (cervical c-1)",
  "Brain - Substantia nigra",
  "Whole Blood"
)

# ============================================================================
# 单个组织的分析函数 (用于并行)
# ============================================================================
process_tissue <- function(tissueSiteDetail, gwasDF, sentinelSnpDF, gwas_name, 
                           output_dir, egenes_cache_dir, force_rerun, max_retries = 10) {
  
  # 清理组织名称用于文件名
  tissue_clean <- str_replace_all(tissueSiteDetail, " |\\(|\\)", "_")
  out_file <- file.path(
    output_dir,
    paste0(gwas_name, "_", tissue_clean, "_colocResultSig.RDS")
  )
  
  # 检查输出文件是否已存在
  if (file.exists(out_file) && !force_rerun) {
    return(list(
      tissue = tissueSiteDetail,
      status = "skipped",
      message = "Result file already exists",
      n_significant = NA
    ))
  }
  
  # 使用替代方案: 直接下载整个组织的eGenes,然后手动筛选
  traitsAll <- NULL
  retry_count <- 0
  
  # 方法1: 尝试使用 xQTLanalyze_getTraits (原始方法)
  while (is.null(traitsAll) && retry_count < 10) {
    traitsAll <- tryCatch({
      options(timeout = 300)
      xQTLanalyze_getTraits(
        sentinelSnpDF,
        detectRange = 1e6,
        tissueSiteDetail = tissueSiteDetail
      )
    }, error = function(e) {
      retry_count <<- retry_count + 1
      if (retry_count < 3) {
        Sys.sleep(10)
      }
      return(NULL)
    })
  }
  
  # 方法2: 如果方法1失败,使用 xQTLdownload_egene 下载整个组织的eGenes
  if (is.null(traitsAll)) {
    # 检查是否有缓存的eGenes数据
    egenes_cache_file <- file.path(egenes_cache_dir, paste0(tissue_clean, "_eGenes.RDS"))
    eGenesAll <- NULL
    
    if (file.exists(egenes_cache_file)) {
      tryCatch({
        eGenesAll <- readRDS(egenes_cache_file)
      }, error = function(e) {
        eGenesAll <<- NULL
      })
    }
    
    # 如果没有缓存或缓存加载失败,则下载
    if (is.null(eGenesAll)) {
      retry_count <- 0
      
      while (is.null(eGenesAll) && retry_count < max_retries) {
        eGenesAll <- tryCatch({
          options(timeout = 300)
          xQTLdownload_egene(
            tissueSiteDetail = tissueSiteDetail
          )
        }, error = function(e) {
          retry_count <<- retry_count + 1
          if (retry_count < max_retries) {
            Sys.sleep(10)
          }
          return(NULL)
        })
      }
      
      # 保存下载的eGenes数据到缓存
      if (!is.null(eGenesAll) && nrow(eGenesAll) > 0) {
        tryCatch({
          saveRDS(eGenesAll, file = egenes_cache_file)
        }, error = function(e) NULL)
      }
    }
    
    if (is.null(eGenesAll) || nrow(eGenesAll) == 0) {
      return(list(
        tissue = tissueSiteDetail,
        status = "failed",
        message = "No eGenes found",
        n_significant = 0
      ))
    }
    
    # 手动筛选: 找到在sentinel SNP附近的基因
    traitsAll <- data.table()
    detectRange <- 1e6
    
    for (i in 1:nrow(sentinelSnpDF)) {
      sentinel_chr <- sentinelSnpDF$chrom[i]
      sentinel_pos <- sentinelSnpDF$position[i]
      sentinel_snp <- sentinelSnpDF$rsid[i]
      
      # 查询该染色体上该位置附近的基因
      genesNearby <- tryCatch({
        xQTLquery_gene(
          chrom = sentinel_chr,
          start = max(1, sentinel_pos - detectRange),
          end = sentinel_pos + detectRange
        )
      }, error = function(e) NULL)
      
      if (!is.null(genesNearby) && nrow(genesNearby) > 0) {
        # 找到这些基因中哪些是eGenes
        eGenesNearby <- eGenesAll[eGenesAll$gencodeId %in% genesNearby$gencodeId, ]
        
        if (nrow(eGenesNearby) > 0) {
          # 构建traitsAll格式的数据
          traits_i <- data.table(
            rsid = sentinel_snp,
            chrom = sentinel_chr,
            position = sentinel_pos,
            gencodeId = eGenesNearby$gencodeId,
            geneSymbol = eGenesNearby$geneSymbol
          )
          traitsAll <- rbind(traitsAll, traits_i)
        }
      }
    }
    
    if (nrow(traitsAll) == 0) {
      return(list(
        tissue = tissueSiteDetail,
        status = "failed",
        message = "No eGenes found near sentinel SNPs",
        n_significant = 0
      ))
    }
  }
  
  # 查询基因信息(添加重试机制)
  genesAll <- NULL
  retry_count <- 0
  
  while (is.null(genesAll) && retry_count < max_retries) {
    genesAll <- tryCatch({
      options(timeout = 300)
      xQTLquery_gene(unique(traitsAll$gencodeId))
    }, error = function(e) {
      retry_count <<- retry_count + 1
      if (retry_count < max_retries) {
        Sys.sleep(10)
      }
      return(NULL)
    })
  }
  
  if (is.null(genesAll) || nrow(genesAll) == 0) {
    return(list(
      tissue = tissueSiteDetail,
      status = "failed",
      message = "No genes found",
      n_significant = 0
    ))
  }
  
  # 初始化结果存储
  colocResultAll <- data.table()
  
  # 遍历基因
  for (i in 1:nrow(genesAll)) {
    
    gene_id <- genesAll[i]$gencodeId
    gene_symbol <- genesAll[i]$geneSymbol
    
    # 下载eQTL数据(添加重试机制,尝试两个数据源)
    eQTL_i <- NULL
    retry_count <- 0
    data_sources <- c("liLab", "eQTL_catalogue")
    
    while (is.null(eQTL_i) && retry_count < max_retries) {
      # 交替使用两个数据源
      current_source <- data_sources[(retry_count %% 2) + 1]
      
      eQTL_i <- tryCatch({
        options(timeout = 300)
        xQTLdownload_eqtlAllAsso(
          gene_id,
          geneType = "gencodeId",
          tissueLabel = tissueSiteDetail,
          withB37VariantId = FALSE,
          data_source = current_source
        )
      }, error = function(e) {
        retry_count <<- retry_count + 1
        if (retry_count < max_retries) {
          Sys.sleep(3)
        }
        return(NULL)
      })
    }
    
    if (is.null(eQTL_i) || nrow(eQTL_i) == 0) {
      next
    }
    
    # 匹配GWAS和eQTL数据
    gwasDF_i <- gwasDF[rsid %in% eQTL_i$rsid]
    if (nrow(gwasDF_i) == 0) {
      next
    }
    
    # 合并位置信息
    eQTL_i <- merge(eQTL_i, gwasDF_i[, .(rsid, chrom, position)], by = "rsid")
    eQTL_i <- eQTL_i[, .(rsid, chrom, position, pValue, maf, beta, se)]
    
    # 清理P值
    gwasDF_i <- clean_pval(gwasDF_i, "pValue")
    eQTL_i <- clean_pval(eQTL_i, "pValue")
    
    # 运行共定位分析
    colocResult_i <- tryCatch({
      xQTLanalyze_coloc_diy(
        gwasDF = gwasDF_i,
        qtlDF = eQTL_i,
        method = "Both"
      )
    }, error = function(e) NULL)
    
    if (!is.null(colocResult_i) && !is.null(colocResult_i$coloc_Out_summary)) {
      colocResult_i <- colocResult_i$coloc_Out_summary
      colocResult_i <- cbind(
        genesAll[i, c("geneSymbol", "gencodeId")],
        colocResult_i
      )
      colocResultAll <- rbind(colocResultAll, colocResult_i)
    }
  }
  
  # 过滤显著结果 (PP.H4.abf > 0.75)
  n_significant <- 0
  if (nrow(colocResultAll) > 0) {
    colocResultsig <- colocResultAll[PP.H4.abf>0.75 & hypr_posterior>0.5][order(-PP.H4.abf)]
    n_significant <- nrow(colocResultsig)
    
    # 保存结果
    saveRDS(colocResultsig, file = out_file)
    
    # 可视化(如果有显著结果) - 在并行中跳过可视化以避免问题
    # 可视化部分可以在主进程中单独处理
  }
  
  return(list(
    tissue = tissueSiteDetail,
    status = "completed",
    message = paste0("Found ", n_significant, " significant colocalization(s)"),
    n_significant = n_significant,
    out_file = out_file
  ))
}

# ============================================================================
# 主分析流程
# ============================================================================
message(paste0(rep("=", 70), collapse = ""))
message("Starting Colocalization Analysis (Parallel Version)")
message("Input directory: ", opt$input_dir)
message("Output directory: ", opt$output_dir)
message("P-value threshold: ", opt$pvalue)
message("Force rerun: ", opt$force)
message("Number of cores: ", opt$ncores)
message(paste0(rep("=", 70), collapse = ""))

# 设置并行环境
n_cores <- min(opt$ncores, detectCores() - 1, length(tissues))
message("Using ", n_cores, " cores for parallel processing")

# 查找所有GWAS文件
gwas_files <- list.files(opt$input_dir, pattern = "smr\\.txt\\.gz$", full.names = TRUE)

if (length(gwas_files) == 0) {
  stop("No files matching pattern '*.smr.txt.gz' found in ", opt$input_dir)
}

message("Found ", length(gwas_files), " GWAS file(s) to process")

# 遍历GWAS文件
for (gwas_path in gwas_files) {
  
  gwas_name <- str_replace(basename(gwas_path), "\\.smr\\.txt\\.gz$", "")
  message("\n", paste0(rep("=", 70), collapse = ""))
  message("Processing GWAS: ", gwas_name)
  message(paste0(rep("=", 70), collapse = ""))
  
  # 读取GWAS文件
  gwasDF <- NULL
  tryCatch({
    gwasDF <- fread(gwas_path, 
                    select = c("SNP", "A1", "A2", "freq", "b", "se", "P", "N", "pos", "chr"))
    setDT(gwasDF)
    message("Loaded ", nrow(gwasDF), " variants from ", basename(gwas_path))
    
  }, error = function(e) {
    message("ERROR: Could not read ", gwas_path, " - ", e$message)
  })
  
  if (is.null(gwasDF) || nrow(gwasDF) == 0) {
    message("Skipping empty or invalid file")
    next
  }
  
  # 过滤只保留rs开头的SNP
  gwasDF <- gwasDF[str_detect(SNP, "^rs")]
  message("After filtering rs SNPs: ", nrow(gwasDF), " variants")
  
  # 检测基因组版本
  version_info <- detect_genome_version(gwasDF)
  
  # 准备数据格式
  gwasDF <- gwasDF[, .(
    rsid = SNP,
    chrom = as.character(chr),
    position = as.integer(pos),
    pValue = as.numeric(P),
    AF = as.numeric(freq),
    beta = as.numeric(b),
    se = as.numeric(se)
  )]
  
  # 清理P值
  gwasDF <- clean_pval(gwasDF, "pValue")
  setindex(gwasDF, rsid)
  
  message("Cleaned data: ", nrow(gwasDF), " variants ready for analysis")
  
  # 获取sentinel SNP
  message("Identifying sentinel SNPs...")
  sentinelSnpDF <- tryCatch({
    xQTLanalyze_getSentinelSnp(
      gwasDF,
      pValueThreshold = opt$pvalue,
      centerRange = 1e6,
      genomeVersion = version_info$version,
      grch37To38 = version_info$grch37To38
    )
  }, error = function(e) {
    message("ERROR in getSentinelSnp: ", e$message)
    return(NULL)
  })
  
  if (is.null(sentinelSnpDF) || nrow(sentinelSnpDF) == 0) {
    message("No sentinel SNPs found, skipping this GWAS file")
    next
  }
  
  message("Found ", nrow(sentinelSnpDF), " sentinel SNP(s)")
  
  # 创建eGenes缓存目录
  egenes_cache_dir <- file.path(opt$output_dir, "egenes_cache")
  if (!dir.exists(egenes_cache_dir)) {
    dir.create(egenes_cache_dir, recursive = TRUE)
  }
  
  # ========================================================================
  # 并行处理组织
  # ========================================================================
  message("\nStarting parallel processing of ", length(tissues), " tissues...")
  
  # 注册并行后端
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 导出必要的函数和变量到worker节点
  clusterExport(cl, c("clean_pval", "process_tissue"), envir = environment())
  
  # 在每个worker上加载必要的包
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(data.table)
      library(stringr)
      library(xQTLbiolinks)
      library(coloc)
      library(hyprcoloc)
    })
  })
  
  # 并行执行
  results <- foreach(
    tissueSiteDetail = tissues,
    .combine = rbind,
    .packages = c("data.table", "stringr", "xQTLbiolinks", "coloc", "hyprcoloc"),
    .errorhandling = "pass"
  ) %dopar% {
    result <- tryCatch({
      process_tissue(
        tissueSiteDetail = tissueSiteDetail,
        gwasDF = gwasDF,
        sentinelSnpDF = sentinelSnpDF,
        gwas_name = gwas_name,
        output_dir = opt$output_dir,
        egenes_cache_dir = egenes_cache_dir,
        force_rerun = opt$force
      )
    }, error = function(e) {
      list(
        tissue = tissueSiteDetail,
        status = "error",
        message = as.character(e),
        n_significant = NA
      )
    })
    
    # 转换为data.frame以便combine
    as.data.frame(result, stringsAsFactors = FALSE)
  }
  
  # 停止并行集群
  stopCluster(cl)
  
  # 汇总结果
  message("\n", paste0(rep("-", 50), collapse = ""))
  message("Summary for GWAS: ", gwas_name)
  message(paste0(rep("-", 50), collapse = ""))
  
  if (is.data.frame(results)) {
    for (i in 1:nrow(results)) {
      status_symbol <- switch(
        as.character(results$status[i]),
        "completed" = "✓",
        "skipped" = "○",
        "failed" = "✗",
        "error" = "!"
      )
      message(sprintf("  %s %s: %s", 
                      status_symbol, 
                      results$tissue[i], 
                      results$message[i]))
    }
    
    # 统计
    n_completed <- sum(results$status == "completed")
    n_skipped <- sum(results$status == "skipped")
    n_failed <- sum(results$status == "failed")
    n_error <- sum(results$status == "error")
    total_significant <- sum(as.numeric(results$n_significant), na.rm = TRUE)
    
    message("\nStatistics:")
    message("  Completed: ", n_completed)
    message("  Skipped (already exists): ", n_skipped)
    message("  Failed: ", n_failed)
    message("  Errors: ", n_error)
    message("  Total significant colocalizations: ", total_significant)
  }
  
  # ========================================================================
  # 后处理：可视化（可选，在主进程中执行）
  # ========================================================================
  message("\nGenerating visualizations for significant results...")
  
  for (tissueSiteDetail in tissues) {
    tissue_clean <- str_replace_all(tissueSiteDetail, " |\\(|\\)", "_")
    out_file <- file.path(
      opt$output_dir,
      paste0(gwas_name, "_", tissue_clean, "_colocResultSig.RDS")
    )
    
    if (file.exists(out_file)) {
      colocResultsig <- readRDS(out_file)
      
      if (nrow(colocResultsig) > 0) {
        outGenes <- tryCatch({
          xQTLquery_gene(colocResultsig$gencodeId)
        }, error = function(e) NULL)
        
        if (!is.null(outGenes) && nrow(outGenes) > 0) {
          outGenes <- merge(
            colocResultsig[, .(gencodeId, PP.H4.abf, candidate_snp, SNP.PP.H4)],
            outGenes[, .(geneSymbol, gencodeId, entrezGeneId, geneType)],
            by = "gencodeId",
            sort = FALSE
          )
          outGenes <- outGenes[geneType == "protein coding"]
          
          if (nrow(outGenes) > 0) {
            tryCatch({
              xQTLvisual_genesExp(
                outGenes$geneSymbol,
                tissueSiteDetail = tissueSiteDetail
              )
            }, error = function(e) {
              message("  Warning: Could not create visualization for ", tissueSiteDetail, " - ", e$message)
            })
          }
        }
      }
    }
  }
  
} # end GWAS loop

message("\n", paste0(rep("=", 70), collapse = ""))
message("Analysis complete!")
message(paste0(rep("=", 70), collapse = ""))