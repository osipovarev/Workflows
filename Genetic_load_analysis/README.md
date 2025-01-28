# Genetic Load Analysis Workflow

---
### Some useful resources
- Great review on genetic load
Dussex et al, 2023: [Purging and accumulation of genetic load in conservation](https://doi.org/10.1016/j.tree.2023.05.008)
- Review: defentitions and how to look for genetic load
Bertorelle et al, 2022: [Genetic load: genomic estimates and applications in non-model animals](https://www.nature.com/articles/s41576-022-00448-x)


This workflow outlines the general steps to perform genetic load analysis. 
There are several ways to identify variant effect on fitness. Conceptually, they can be devided into two main classes: 
1) based on gene annotation;
2) based on conservation.\
The examples of **Class 1** tools include  **snpEff** and **ANNOVAR**.\
An example of **Class 2** tool is **GERP**.

Below I describe how to use all three tools to annotate variant effect.



---

## 1. Variant Annotation with snpEff

- snpEff output annotations explained:
[table describing types of variants](https://pcingola.github.io/SnpEff/snpeff/inputoutput/#eff-field-vcf-output-files)


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

Publicaiton:
[ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC2938201/)

- link to [download ANNOVAR](http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz)

- [Code](http://www.openbioinformatics.org/annovar/)

- [tutorial for non-human species](https://annovar.openbioinformatics.org/en/latest/user-guide/gene/#create-your-own-gene-definition-databases-for-non-human-species)


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



---

## 3. Genetic Constraint Scoring with GERP

GERP = Genomic Evolutionary Rate Profiling
Publication: Davydov, ... Sidow, 2010: [Identifying a High Fraction of the Human Genome to be under Selective Constraint Using GERP++](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001025)

### Main idea:
if a site is under purifying selection, it will be conserved across large evolutionary distances.
We can use GERP conservation scores to annotate variant effect: higher conservation of a genomic position means higher impact of a variant in that position.

### What is GERP score?
The GERP score is defined as the reduction in the number of substitutions in the multi-species sequence alignment compared to the neutral expectation.
For example, a GERP score of 4 would mean there are 4 fewer substitutions at a particular site than what is expected based on the neutral rate of evolution across the phylogeny.


### Input required for GERP
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


---

****
## Other methods to predict consequence of sequence variation

**Methods that rely on conservation only:**
- SIFT (https://academic.oup.com/nar/article/31/13/3812/2904131) (requires MSA!)
- PROVEAN (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0046688)
- EVmutation (https://www.nature.com/articles/nbt.3769)
- phyloP (https://academic.oup.com/bib/article/12/1/41/244593?login=false)



**Methods that rely on annotation+ conservation:**
- LIST: [Improved measures for evolutionary conservation that exploit taxonomy distances](https://www.nature.com/articles/s41467-019-09583-2)
- PolyPhen-2 (https://doi.org/10.1002/0471142905.hg0720s76) (requires MSA!)
- CADD (https://www.nature.com/articles/ng.2892)
- Eigen (https://www.nature.com/articles/ng.3477)
- DANN (https://doi.org/10.1093/bioinformatics/btu703)
- fitCons (https://www.nature.com/articles/ng.3196)