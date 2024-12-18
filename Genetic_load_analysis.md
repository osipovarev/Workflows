
# 1. Genetic load analysis

### Count total number of variants by impact
```
WDIR=

for f in $(ls *.g.genes.txt); \
do i=${f%.g.genes.txt}; \
	high=$(cut -f5 $f | tail -n +3 | sum_stdin.py); \
	low=$(cut -f6 $f | tail -n +3 | sum_stdin.py); \
	mod=$(cut -f7 $f | tail -n +3 | sum_stdin.py); \
	modif=$(cut -f5,6,7,8 $f | awk '$0=$4-$1-$2-$3{print}' | tail -n +3 | sum_stdin.py); \
	 echo -e "$i\t$high\t$low\t$mod\t$modif"; \
done > $WDIR/all_samples.total_variant_by_impact.tsv

```

### Count number of homozygotes and heterozygotes by impact
```
WDIR=
SAMPLES=

### wrote jobs.bcftools_by_impact.sh
sbatch --array=1-153 jobs.bcftools_by_impact.sh


for SAMPLE in $(cat $SAMPLES); \
do \
	for STATE in hom het; \
	do \
		for IMPACT in HIGH MODERATE LOW; \
		do \
			NUMBER=$(grep -v ^# $WDIR/$IMPACT.$STATE.$SAMPLE.snpEff.vcf | wc -l); \
			echo -e "$SAMPLE\t$IMPACT\t$STATE\t$NUMBER"; \
		done; \
	done; \
done > all_samples.hom_het_by_impact.tsv

### add LOF (subset of HIGH) variants separately
IMPACT=HIGH

for SAMPLE in $(cat $SAMPLES); \
do \
        for STATE in hom het; \
        do \
                NUMBER=$(grep -v ^# $WDIR/$IMPACT.$STATE.$SAMPLE.snpEff.vcf | grep -w LOF | wc -l); \
        echo -e "$SAMPLE\tLOF\t$STATE\t$NUMBER"; \
        done; \
done >> all_samples.hom_het_by_impact.tsv

```


### Re-count number of homozygotes splitting into hom reference (0/0) and hom alternate (1/1)
```

awk '{OFS = "\t"; $6=$5-$4; print}' all_samples.hom_het_by_impact_with_ref.tsv | awk '$6!=0{print}' | cut -f1-3,6 | sed 's/hom/hom_ref/' > hom_ref.by_impact.tsv

cut -f1-4 all_samples.hom_het_by_impact_with_ref.tsv | sed 's/hom/hom_alt/' > hom_alt.by_impact.tsv

cat hom_ref.by_impact.tsv hom_alt.by_impact.tsv > all_samples.hom_het_by_impact.ref_alt.tsv
```


# 2. Inbreeding and genetic load

### Calculation of F_ROH for each individual of leatherback turtle

genome size = 2,164,762,090 bp

```
WDIR=
SAMPLES=

sbatch --array=1-153 jobs.bcftools_roh.sh


for SAMPLE in $(cat $SAMPLES); \
do \
	L_ROH=$(grep ^RG ROH_per_individual/$SAMPLE.roh.txt | awk '{sum += $6} END {print sum}'); \
	echo -e "$SAMPLE\t$L_ROH"; \
done > all_samples.L_ROH.tsv
```

