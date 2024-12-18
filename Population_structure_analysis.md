# Population structure analysis

This recipe descibes how to run population structure analysis with PLINK and ADMIXTURE

helpful resources:
1) [Exploring population structure with admixture models and principal components analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8722024/)


## 0. Load modules
```
module load uri/main
module load PLINK/2.00a3.7-gfbf-2023a # loads plink v1.9

module load bcftools/1.19

```

## 1. PLINK

### 1.1. Convert to bed; prepare input
```
VCF=
PREFIX=

plink --vcf $VCF --make-bed --double-id --allow-extra-chr --chr-set 95 no-xy --threads 12 --out $PREFIX
```
This takes ~40 mins for 10M variants
NB: argument `--chr-set 95 no-xy` is needed when chromosomes don't have numberic names (1, 2, 3 ..), but for exmaple SUPER_01


### Optional: filter empty positions
Sometimes PLINK produces input with positions where all values are missing;
If this happens, you can filter them out running:
```
plink --bfile $PREFIX --geno 0.99 --maf 0.001 --chr-set 95 no-xy --make-bed --out filt.$PREFIX
```


### 1.2. Filter out SNPs to remove linkage disequilibrium (LD)
SNPs in high LD with each other contain redundant information and can disproportionately influence the results of the population structure analysis. A standard approach to address this issue is to filter out SNPs based on pairwise LD to produce a reduced set of more independent markers.

```
# Remove linked SNPs (LD pruning) with r2 > 0.1:
plink --bfile $PREFIX --indep-pairwise 50 10 0.1 --chr-set 95 no-xy --make-bed --threads 32 --out ld_pruned_0.1.$PREFIX

# Remove linked SNPs (LD pruning) with r2 > 0.2:
plink --bfile $PREFIX --indep-pairwise 50 10 0.2 --chr-set 95 no-xy --make-bed --threads 32 --out ld_pruned_0.2.$PREFIX
```



## 2. PCA

### 2.A. PCA from PLINK
```
plink --bfile turtles_ld_pruned --allow-extra-chr --chr-set 95 no-xy --threads 32  --pca --out pca.$PREFIX
```
NB: this fails with Illegal instruction (core dumped). :(


### 2.B. PCA from distance matrix
```
plink --bfile turtles_ld_pruned --allow-extra-chr --chr-set 95 no-xy --threads 32  --distance-matrix --out pca.$PREFIX
```



## 3. Admixture
To run ADMIXTURE, we need to give the program only two things: 
1) the SNP data (.bed file created with PLINK) and
2) a number of ancestral populations (K) to estimate ancestry proportions for.

### 3.1. Make sbatch script:
```
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
PREFIX=
admixture --cv ld_pruned_0.2.$PREFIX.bed $K
```

### 3.2. Run ADMIXTURE 
This runs ADMIXTURE for ancestral populations (K) from 2 to 8.
```
sbatch --array=2-8 jobs.admixture.sh
```
NB: To assess what the best value of K is, run ADMIXTURE with cross-validation by adding the --cv flag to the command.

Check the logs to see which K value corresponds to the lowest CV error.



## Alternative: fastmixture

You can check out the newer afster version of ADMIXTURE - fastmixture:
Santander, et al 2024: [Faster model-based estimation of ancestry proportions](https://www.biorxiv.org/content/10.1101/2024.07.08.602454v3)

