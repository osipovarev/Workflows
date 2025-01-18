# Genetic Load Analysis Workflow

This workflow outlines the general steps to perform genetic load analysis. 
There are several ways to identify variant effect on fitness. Conceptually, they can be devided into two main classes: 
1) based on gene annotationl and
2) based on conservation.
The examples of **Class 1** tools include  **snpEff** and **ANNOVAR**.
An example of **Class 2** tool is **GERP**.

Below I describe how to use all three tools to annotate variant effect.

---

## 1. Input Preparation

### Step 1.1: Prepare Required Data
- **Genome Reference**: Obtain the genome sequence (FASTA) and annotation file (GTF).
- **Variants**: Ensure your variant data is in VCF format.

### Step 1.2: Organize Working Directory
Structure your working directory for the tools:
- A folder for the genome and annotation files.
- Subdirectories for outputs from snpEff, ANNOVAR, and GERP.

---

## 1. Variant Annotation with snpEff


### 1.1. Prepare input: VCF, genome, annotation

1. In `data/` directory, put genomic fasta file and annotation.gtf file (or links to these fiels)

the structure that snpEff understands is:
DB_NAME=[rDerCor1.pri.cur.20210524] - example name

data/genomes/$DB_NAME.fa
data/$DB_NAME/genes.gtf

2. Add a line corresponding to genome and annotation to the snpEff.config file, e.g.:
```
$DB_NAME.genome : $DB_NAME
```

3. Prepare combined VCF with all samples you want to annotate.



### 1.2. Annotate Variants with snpEff

1. Use snpEff to classify variants by impact (HIGH, MODERATE, LOW, MODIFIER).
(Edit job.snpEff_combined.sh providing your input data)

  ```bash
  java -jar snpEff.jar build -gff3 genome
  java -jar snpEff.jar genome input.vcf > annotated.vcf
  ```

2. Split combined annotated VCF by impact; then split by sample
(Edit job.by_impact_conbined.sh and jobs.bcftools_split_samples.sh providing your input data)

  ```bash
sbatch job.by_impact_conbined.sh

sbatch --array=1-153 jobs.bcftools_split_samples.sh HIGH
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODERATE
sbatch --array=1-153 jobs.bcftools_split_samples.sh LOW
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODIFIER
  ```

### 1.3. Count number of snps and indels of homozygotes and heterozygotes by impact
(Edit count_var_by_impact.sh)
```
bash count_var_by_impact.sh
```

### 1.4. Make combined table with total number of homozygotes
```
paste <(grep -v -w LOF all_samples.hom_het_by_impact.tsv ) <(cut -f4 snpEff_combined/all_samples.hom_het_by_impact.tsv)  > all_samples.hom_het_by_impact_with_ref.tsv
```



## 2. Annotate Variants with ANNOVAR


### 2.1. Prepare input: genome, annotation, VCF

Prepare a custom database from the genome annotation.
This custom database should be in the annovar/ repository that you downloaded.

DB_NAME=[turtledb]
put genes.txt in annovar/$DB_NAME


### 2.2. Annotate variants
(Edit job.annovar_conbined.sh adding your input data)
  ```bash
  sbatch job.annovar_conbined.sh
  ```

Since effect annotations by snpEff and ANNOVAR differ, here is the table of correspondance between ANNOVAR and snpEff: 

| ANNOVAR annotation | snpEff annotation |
|--------------------|-------------------|
| frameshift_insertion          |   HIGH | 
| frameshift_deletion          | HIGH |
| frameshift_block_substitution |   HIGH |
| stopgain                     |  HIGH |
| stoploss                     |  HIGH |
| nonframeshift_insertion          | MODERATE |
| nonframeshift_deletion          |  MODERATE |
| nonframeshift_block_substitution  |   MODERATE |
| nonsynonymous_SNV           | MODERATE |
| synonymous_SNV                  |  LOW |
| unknown          | MODIFIER |

Therefore, the next step explains how to re-annotate ANNOVAR effects to match snpEff annotation

### 2.3. Re-Classify and split variants by impact
  ```bash
  sbatch job.by_impact_annovar.sh

## fix header
VCF=...

echo '##INFO=<ID=.,Number=0,Type=Flag,Description="Placeholder for undefined INFO">' > header.txt

for impact in HIGH MODERATE LOW MODIFIER; \
do \
    bcftools annotate -h header.txt $impact.$VCF -O z -o fixed.vcf.gz; \
    mv fixed.vcf.gz $impact.$VCF; \
done
```


### 2.4. Split by sample and by state (hom/het)
(Edit jobs.bcftools_split_samples.sh)
```
sbatch --array=1-153 jobs.bcftools_split_samples.sh HIGH
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODERATE
sbatch --array=1-153 jobs.bcftools_split_samples.sh LOW
sbatch --array=1-153 jobs.bcftools_split_samples.sh MODIFIER
```

### 2.5. Count number of snps and indels of homozygotes and heterozygotes by impact
(Edit count_var_by_impact.sh)
```
bash count_var_by_impact.sh
```





## 3. Genetic Constraint Scoring with GERP

GERP uses conservation to annotate variant effect: higher conservation of a genomic position means higher impact of a variant in that position.

To identify conserved position, GERP needs a multiple sequence alignment across species with large evoltuionary distance. This recipe does not describe how to make this alignemnt.

To install gerp, run
```
mamba install bioconda::gerp
```


### 3.1. Prepare input: MAF, annotation

To make GERP run faster, split multiple sequence alignemnt in MAF format into chromosomes/scaffolds.


### 3.2: Run GERP Analysis

jobs.filter_maf.sh filters MAF alignemnt

jobs.gerp.sh does two things:
1) it calculates GERP score for each position in the MAF alignment;
2) it merges stretches of high-scoring position into 'conserved elements'

GERP tends to be very greedy with how it merges position into conserved elements, so if you are looking to use the elements, consider breaking them up using information from other tools, like PhastCons or PhyloP.

  ```bash
  sbatch --array=1-40 jobs.filter_maf.sh
  sbatch --array=1-40 jobs.gerp.sh
  ```

### 3.3. Get GERP summary stats in histograms
  ```bash
# distribution of scores
cat mafs/*rates     | cut -f2 | python3 hist_stdin.py 20 gerp.rates.hist.pdf

# distribution of sum scores per conserved element
cat mafs/*maf.rates.elems | cut -f4 | python3 hist_stdin.py 100 elem.rates.hist.pdf

# distribution of lenghts of conserved elements
cat mafs/*maf.rates.elems | cut -f6 | python3 hist_stdin.py 50 elem.length.hist.pdf
  ```


### 3.4. Separate GERP results into coding and non-coding regions
- Use bedtools to extract coding and non-coding regions:
  ```bash
  sbatch --array=1-40 jobs.bedtools.sh
  ```
- Check distribution of scores in coding and non-coding regions:
```
cat mafs/coding*bed | cut -f4 | python3 hist_stdin.py 50 coding.gerp.rates.hist.pdf

cat mafs/noncoding*bed | cut -f4 | python3 hist_stdin.py 50 noncoding.gerp.rates.hist.pdf
```


### 3.5. Get positions with high GERP score (conserved). Overlap with annotation
```
ANNO=...
EXONS=...

## Split annotation into coding exons
bed12ToBed6 -i $ANNO > $EXONS


## Extract positions with GERP score >= 1.0
### coding
bedtools merge -i <(cat coding.*.maf.rates.bed | awk '$4>=1{print }' | sort -k1,1 -k2,2n) > ../coding.derCor.gerp_1.0.bed

### noncoding
bedtools merge -i <(cat noncoding.*.maf.rates.bed | awk '$4>=1{print }' | sort -k1,1 -k2,2n) > ../noncoding.derCor.gerp_1.0.bed
```

