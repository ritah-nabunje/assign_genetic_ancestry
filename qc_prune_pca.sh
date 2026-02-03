#!/bin/bash

# add job scheduler details here if needed

#-----------------
# 1. qc_prune_pca.sh
# Perform QC on reference and study genotype data, harmonise variants,
#   select shared LD-pruned SNPs, compute PCA in a reference panel,
#   and project study samples into the same PC space.
# Requires:
#   - PLINK 2.0
#   - BED/BIM/FAM for both ref and study data
#   - High-LD region file in PLINK's --exclude range format
# Output:
# 	- reference PCs (eigenvec and eigen values)
#   - study PCs 
#----------------------------------------------

WORKING_DIR="/path/to/dir"

# reference genotype data (PLINK prefix)
REF="/path/to/reference"
REF_NAME="ref"

# study genotype data (PLINK prefix)
STUDY="/path/to/genotypes"

# regions of extended LD to exclude (e.g. HLA)
HIGH_LD_REGIONS="/path/to/exclude.txt"

cd "${WORKING_DIR}"
module load plink2

# --------------QC on ref
plink2 \
	--bfile "${REF}" \
	--allow-no-sex \
	--autosome \
	--geno 0.05 \
	--mind 0.05 \
	--maf 0.01 \
	--hwe 1e-6 \
	--exclude range "${HIGH_LD_REGIONS}" \
	--make-bed \
	--out "${REF_NAME}.qc"

# remove strand-ambiguous SNPs (A/T, C/G) from ref
awk 'BEGIN {OFS="\t"} ($5$6 == "GC" || $5$6 == "CG" || $5$6 == "AT" || $5$6 == "TA") {print $2}' "${REF_NAME}.qc.bim" > "${REF_NAME}.atcg_snps"

plink2 \
	--bfile "${REF_NAME}.qc" \
	--allow-no-sex \
	--exclude "${REF_NAME}.atcg_snps" \
	--make-bed \
	--out "${REF_NAME}.qc.no_atcg"

# check if variant IDs are same in ref and study data, if not, 
# harmonise variant IDs in study data
# Format: chr:pos:ref:alt = @:#:\$r:\$a 
# match the reference variant ID convention
plink2 \
	--bfile "${STUDY}" \
	--allow-no-sex \
	--set-all-var-ids @:#:\$r:\$a \
	--new-id-max-allele-len 50 missing \
	--make-bed \
	--out "${STUDY}.ids"

#----------- QC on study data
plink2 \
	--bfile "${STUDY}.ids" \
	--allow-no-sex \
	--geno 0.05 \
	--autosome \
	--maf 0.01 \
	--hwe 1e-6 \
	--mind 0.05 \
	--exclude range "${HIGH_LD_REGIONS}" \
	--make-bed \
	--out "${STUDY}.qc"

# identify shared variants
awk '{print $2}' "${STUDY}.qc.bim" > study.snps
awk '{print $2}' "${REF_NAME}.qc.no_atcg.bim" > ref.snps

sort study.snps ref.snps | uniq -c | awk '$1==2 {print $2}' > shared.variants

# extract shared variants in ref
plink2 \
	--bfile "${REF_NAME}.qc.no_atcg" \
	--allow-no-sex \
	--extract shared.variants \
	--make-bed \
	--out "${REF_NAME}.shared"

# keep independent variants 
# LD pruning on reference 
plink2 \
	--bfile "${REF_NAME}.shared" \
	--allow-no-sex \
	--indep-pairwise 50 5 0.3 \
	--out "${REF_NAME}.prune"

# extract ref prune.in variants from ref
plink2 \
	--bfile "${REF_NAME}.shared" \
	--allow-no-sex \
	--extract "${REF_NAME}.prune.prune.in" \
	--make-bed \
	--out "${REF_NAME}.pruned"

# Extract ref prune.in variants from study data
plink2 \
	--bfile "${STUDY}.qc" \
	--allow-no-sex \
	--extract "${REF_NAME}.prune.prune.in" \
	--make-bed --out "${STUDY}.pruned"

# PCA in reference and projection of study samples
# See: https://www.cog-genomics.org/plink/2.0/score#pca_project
# could do this differently though - combining the ref and study data into one set, then calcualte PCs
# you'd have to separate the ref and study PCs for the next steps

## calculate ref PCs and PC loadings 
plink2 \
	--bfile "${REF_NAME}.pruned" \
	--allow-no-sex \
	--freq counts \
	--pca allele-wts \
	--out "${REF_NAME}.pca"

# project study samples into the reference PC space
plink2 \
	--bfile "${STUDY}.pruned" \
	--allow-no-sex \
	--read-freq "${REF_NAME}.pca.acount" \
	--score "${REF_NAME}.pca.eigenvec.allele" 2 5 header-read no-mean-imputation \
	variance-standardize \
	list-variants \
	--score-col-nums 6-15 \
	--out "${STUDY}.pca_projected"

# it may help to visualise the reference and study PCA together first before you continue. 
