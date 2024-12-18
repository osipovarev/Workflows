# Population Structure Analysis

This recipe describes how to run population structure analysis using PLINK and ADMIXTURE.

## Helpful Resources
- [Exploring population structure with admixture models and principal components analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8722024/)



## 0. Load Required Modules

```bash
module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a  # loads plink v1.9
module load bcftools/1.19
```


## 1. PLINK Processing

### 1.1. Convert VCF to PLINK Binary Format and Prepare Input

```bash
VCF=  # Set your VCF file path
PREFIX=  # Set your output prefix

plink --vcf $VCF \
      --make-bed \
      --double-id \
      --allow-extra-chr \
      --chr-set 95 no-xy \
      --threads 12 \
      --out $PREFIX
```

**Notes:**
- Processing takes approximately 40 minutes for 10M variants
- The `--chr-set 95 no-xy` argument is necessary when chromosomes have non-numeric names (e.g., SUPER_01 instead of 1, 2, 3)

### 1.2. Optional: Filter Empty Positions

If PLINK produces input with positions where all values are missing, filter them out:

```bash
plink --bfile $PREFIX \
      --geno 0.99 \
      --maf 0.001 \
      --chr-set 95 no-xy \
      --make-bed \
      --out filt.$PREFIX
```

### 1.3. Linkage Disequilibrium (LD) Pruning

SNPs in high LD contain redundant information and can disproportionately influence population structure analysis. Use the following commands to create LD-pruned datasets:

```bash
# LD pruning with r2 > 0.1
plink --bfile $PREFIX \
      --indep-pairwise 50 10 0.1 \
      --chr-set 95 no-xy \
      --make-bed \
      --threads 32 \
      --out ld_pruned_0.1.$PREFIX

# LD pruning with r2 > 0.2
plink --bfile $PREFIX \
      --indep-pairwise 50 10 0.2 \
      --chr-set 95 no-xy \
      --make-bed \
      --threads 32 \
      --out ld_pruned_0.2.$PREFIX
```

***
## 2. Principal Component Analysis (PCA)

### 2.A. PCA from PLINK (Potentially Problematic)

```bash
plink --bfile $PREFIX \
      --allow-extra-chr \
      --chr-set 95 no-xy \
      --threads 32 \
      --pca \
      --out pca.$PREFIX

# Note: This may fail with an "Illegal instruction" error
```

### 2.B. PCA from Distance Matrix

```bash
plink --bfile $PREFIX \
      --allow-extra-chr \
      --chr-set 95 no-xy \
      --threads 32 \
      --distance-matrix \
      --out pca.$PREFIX
```


***
## 3. ADMIXTURE Analysis

ADMIXTURE requires two inputs:
1. SNP data in .bed format (created with PLINK)
2. Number of ancestral populations (K) to estimate ancestry proportions

### 3.1. Create SLURM Batch Script (jobs.admixture.sh)

```bash
#!/bin/bash
#SBATCH -J admixture
#SBATCH -o ./logs/admixture/admixture_%a_%A.log
#SBATCH -e ./logs/admixture/admixture_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 24:00:00
#SBATCH -c 12
#SBATCH --threads-per-core=1
#SBATCH --mem-per-cpu=20000
#SBATCH --array=2-8
#SBATCH --mail-type=END

K=$SLURM_ARRAY_TASK_ID
PREFIX=  # Set your prefix
admixture --cv ld_pruned_0.2.$PREFIX.bed $K
```

### 3.2. Run ADMIXTURE

```bash
sbatch --array=2-8 jobs.admixture.sh
```

**Notes:**
- ADMIXTURE runs with cross-validation (--cv flag)
- Check logs to identify the K value with the lowest cross-validation error



## Alternative: FastMixture

Consider using FastMixture, a newer and faster version of ADMIXTURE:
- Reference: Santander, et al. 2024 - [Faster model-based estimation of ancestry proportions](https://www.biorxiv.org/content/10.1101/2024.07.08.602454v3)
