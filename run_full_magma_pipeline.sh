#!/bin/bash
set -e

# --- 路径配置 ---
IN_DIR="/home/data/t070502/R/GWAS/vcf/converted_gwas" # 使用绝对路径，不要用 ~PREP_DIR="./magma_input"
RESULT_DIR="./magma_out_without_MHC"
PREP_DIR='./magma_input_without_MHC'
BIM_REF="./g1000_eur.bim"
GENE_LOC="./NCBI37.3.noMHC.gene.loc"  # 注意：既然是 map 到 hg19，请务必用 37 的 loc
REF_BFILE="./g1000_eur"        # MAGMA 参考面板前缀

mkdir -p ${PREP_DIR} ${RESULT_DIR}

# 这样 Bash 会把 * 展开，循环给 Python 处理
for f in ${IN_DIR}/*.smr.txt.gz; do
    base=$(basename "${f}" .smr.txt.gz)
    
    # 如果最后的 MAGMA 结果已经存在，则整个循环跳过该文件
    if [ -f "${RESULT_DIR}/${base}_magma_res.genes.out" ]; then
        echo "ALREADY DONE: ${base}, skipping..."
        continue
    fi

    # 1. Python 预处理 (Python 内部也有跳过逻辑)
    python3 check_and_prep_no_mhc.py "$f" "${BIM_REF}" "${PREP_DIR}/${base}.ready"

    # 如果 Python 没生成文件（比如因为某些原因失败了），不跑 MAGMA
    [ -f "${PREP_DIR}/${base}.ready" ] || continue

    # 2. MAGMA Annotation
    ./magma --annotate --snp-loc "${PREP_DIR}/${base}.ready" --gene-loc "${GENE_LOC}" --out "${RESULT_DIR}/${base}_annot"

    # 3. MAGMA Gene-test
    ./magma --bfile "${REF_BFILE}" --pval "${PREP_DIR}/${base}.ready" use=SNP,P ncol=N \
          --gene-annot "${RESULT_DIR}/${base}_annot.genes.annot" --out "${RESULT_DIR}/${base}_magma_res"
done