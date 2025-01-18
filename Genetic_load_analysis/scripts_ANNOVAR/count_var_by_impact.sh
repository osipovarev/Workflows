WDIR=/nese/meclab/Katya/ANNOVAR_turtles/
SCRATCH=/scratch/workspace/lkomoroske_umass_edu-Leatherback_WGRv2/
SAMPLES=$SCRATCH/Dc_SNPs_May24_samplelist.txt
OUTFILE=all_samples.snp_indel_hom_het_by_impact.tsv

module load bcftools/1.19

for SAMPLE in $(cat $SAMPLES); 
do 
    for STATE in hom het; 
    do 
	for IMPACT in HIGH MODERATE LOW MODIFIER; 
        do 
            for i in snps indels; 
            do 
		if [ $STATE == "hom" ]; 
		then 
			NUMBER=$(bcftools view -v ${i} -i 'GT=="0/0"' $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.vcf.gz | grep -v ^# |  wc -l);
			echo -e "$SAMPLE\t$IMPACT\t${STATE}_ref\t$NUMBER\t$i";
			NUMBER=$(bcftools view -v ${i} -i 'GT!="0/0"' $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.vcf.gz | grep -v ^# |  wc -l);
		else
                	NUMBER=$(bcftools view -v ${i} $WDIR/split_by_sample/$IMPACT.$STATE.$SAMPLE.vcf.gz | grep -v ^# |  wc -l); 
		fi;
                echo -e "$SAMPLE\t$IMPACT\t$STATE\t$NUMBER\t$i"; \
            done; 
        done; 
    done; 
done > $OUTFILE
